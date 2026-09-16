//! The comparison, drawn as it fills.
//!
//! Every row a run will produce is printed once, empty, and each is rewritten
//! in place as its measurements land. A pipeline's own progress is not shown:
//! what a reader wants from a benchmark is the numbers, and a list of the
//! commands that produced them is the table's job.
//!
//! Falls back to one line per finished row where stdout is not a terminal, so
//! piping the run to a file still leaves something to read.

use std::io::{IsTerminal, Write};
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::{Arc, Mutex};
use std::thread;
use std::time::{Duration, Instant};

use michi::{Item, Sink, Step};

use crate::ms;

// --- colour

const RESET: &str = "\x1b[0m";
const BOLD: &str = "\x1b[1m";
const DIM: &str = "\x1b[2m";
const GREEN: &str = "\x1b[32m";
const YELLOW: &str = "\x1b[33m";
const RED: &str = "\x1b[31m";
const CYAN: &str = "\x1b[36m";
const MAGENTA: &str = "\x1b[35m";

/// One frame per tick of the painter, so something moves even while a four
/// second run is in flight.
const SPINNER: [char; 10] = ['⠋', '⠙', '⠹', '⠸', '⠼', '⠴', '⠦', '⠧', '⠇', '⠏'];

/// How often the board repaints itself between events.
const TICK: Duration = Duration::from_millis(80);

/// One colour per row of a group, so a row is recognisable before it is read.
//
// by position rather than by name: under `cli` every
// operation lists sail, then seqkit, then an esl-* miniapp,
// so this gives each tool the same colour throughout without
// the board having to know any of their names
fn colour_at(n: usize) -> &'static str {
    match n {
        0 => CYAN,
        1 => MAGENTA,
        _ => YELLOW,
    }
}

// --- the rows

/// One tool's row under one operation.
struct Row {
    op: &'static str,
    tool: &'static str,

    /// Every wall time reported so far, one per repeat.
    walls: Vec<f64>,
    sys_s: f64,
    max_rss_kb: i64,
}

impl Row {
    fn stats(&self) -> Option<(f64, f64)> {
        let best = self.walls.iter().copied().fold(f64::INFINITY, f64::min);
        if !best.is_finite() {
            return None;
        }

        let n = self.walls.len();
        let mean = self.walls.iter().sum::<f64>() / n as f64;
        let stddev = if n < 2 {
            0.0
        } else {
            (self.walls.iter().map(|w| (w - mean).powi(2)).sum::<f64>() / (n - 1) as f64).sqrt()
        };

        Some((best, stddev))
    }
}

/// What the pipeline is doing right now, for the line under the board.
struct Running {
    label: String,
    line: String,
    since: Instant,
}

struct State {
    rows: Vec<Row>,

    /// What the ratio column is against, for its heading.
    subject: &'static str,

    /// The command in flight, redrawn every tick so its clock runs.
    running: Option<Running>,

    /// When the whole run started, for the total on the status line.
    started: Instant,

    /// Why each operation's tools are doing the same work, printed under the
    /// board because a reader who cannot check the claim cannot read the rows.
    claims: Vec<(&'static str, &'static str)>,
    repeat: usize,
    header: String,
    tty: bool,

    /// Lines the last paint left on screen, which the next one moves back over.
    painted: usize,
}

/// The board, shared with the pipeline that fills it.
#[derive(Clone)]
pub(crate) struct Board(Arc<Mutex<State>>);

/// The repainting thread, stopped when the run ends.
pub(crate) struct Painter {
    stop: Arc<AtomicBool>,
    handle: Option<thread::JoinHandle<()>>,
}

impl Painter {
    /// Stop repainting and wait for the thread to let go of stdout.
    pub(crate) fn stop(mut self) {
        self.stop.store(true, Ordering::Relaxed);
        if let Some(h) = self.handle.take() {
            let _ = h.join();
        }
    }
}

impl Board {
    /// A row per tool per operation, in the order they will run.
    pub(crate) fn new(
        plan: &[(&'static str, &'static str, Vec<&'static str>)],
        header: String,
        subject: &'static str,
        repeat: usize,
    ) -> Board {
        let mut rows = Vec::new();
        let claims = plan.iter().map(|(op, why, _)| (*op, *why)).collect();

        for (op, _, tools) in plan {
            for tool in tools {
                rows.push(Row {
                    op,
                    tool,
                    walls: Vec::new(),
                    sys_s: 0.0,
                    max_rss_kb: 0,
                });
            }
        }

        Board(Arc::new(Mutex::new(State {
            rows,
            subject,
            running: None,
            started: Instant::now(),
            claims,
            repeat,
            header,
            tty: std::io::stdout().is_terminal(),
            painted: 0,
        })))
    }

    /// Draw the empty board and keep it moving.
    //
    // a thread rather than a repaint per event: the events
    // are seconds apart, and a clock that only ticks when
    // something finishes is not a clock
    pub(crate) fn start(&self) -> Painter {
        let mut state = self.0.lock().unwrap();

        if !state.tty {
            println!("{}", state.header);

            return Painter {
                stop: Arc::new(AtomicBool::new(true)),
                handle: None,
            };
        }

        state.started = Instant::now();
        paint(&mut state);
        drop(state);

        let stop = Arc::new(AtomicBool::new(false));
        let shared = self.0.clone();
        let flag = stop.clone();

        let handle = thread::spawn(move || {
            while !flag.load(Ordering::Relaxed) {
                thread::sleep(TICK);

                if let Ok(mut state) = shared.lock() {
                    paint(&mut state);
                }
            }
        });

        Painter {
            stop,
            handle: Some(handle),
        }
    }

    /// Leave the finished board in the scrollback, with what each row rests on.
    pub(crate) fn finish(&self) {
        let mut state = self.0.lock().unwrap();
        state.running = None;

        if state.tty {
            paint(&mut state);
        }

        println!();
        for (op, why) in &state.claims {
            let line = format!("  {DIM}{op:<11} {why}{RESET}");
            println!("{}", if state.tty { line } else { plain(&line) });
        }
        println!();
    }
}

impl Sink for Board {
    fn item_start(&mut self, _step: &Step, _at: usize, item: Item<'_>) -> anyhow::Result<()> {
        let mut state = self.0.lock().unwrap();

        state.running = Some(Running {
            label: item.label(),
            line: item.line().unwrap_or_default(),
            since: Instant::now(),
        });

        Ok(())
    }

    fn item_done(&mut self, _step: &Step, _at: usize, item: Item<'_>) -> anyhow::Result<()> {
        // the build and the warmup carry no bench field,
        // which is what keeps them off the board
        let (Some(op), Some(bench)) = (item.fields().get("op"), item.fields().get("bench")) else {
            return Ok(());
        };
        let Some(t) = item.status().timing() else {
            return Ok(());
        };

        let mut state = self.0.lock().unwrap();

        let repeat = state.repeat;
        let Some(row) = state
            .rows
            .iter_mut()
            .find(|r| r.op == op && r.tool == bench)
        else {
            return Ok(());
        };

        // the fastest run's kernel time, kept with it: two
        // numbers from two different runs describe nothing
        if row.walls.iter().copied().fold(f64::INFINITY, f64::min) > t.wall_s {
            row.sys_s = t.sys_s.unwrap_or(0.0);
        }
        row.walls.push(t.wall_s);
        row.max_rss_kb = row.max_rss_kb.max(t.max_rss_kb.unwrap_or(0));

        let full = row.walls.len() == repeat;
        let (op, tool) = (row.op, row.tool);

        if state.tty {
            paint(&mut state);
        } else if full {
            println!("{}", plain(&line_for(&state, op, tool)));
        }

        Ok(())
    }
}

// --- drawing

/// The first row of `op`, which every ratio in the group is against.
//
// first rather than fastest: the question a group answers is
// what the others cost relative to the one being measured,
// and that is the one written down first
fn subject_of(state: &State, op: &str) -> Option<f64> {
    state
        .rows
        .iter()
        .find(|r| r.op == op)
        .and_then(|r| r.stats().map(|(w, _)| w))
}

/// Where `tool` sits among `op`'s rows.
fn place_of(state: &State, op: &str, tool: &str) -> usize {
    state
        .rows
        .iter()
        .filter(|r| r.op == op)
        .position(|r| r.tool == tool)
        .unwrap_or(0)
}

fn line_for(state: &State, op: &str, tool: &str) -> String {
    let row = state
        .rows
        .iter()
        .find(|r| r.op == op && r.tool == tool)
        .expect("a row that was just written");

    let place = place_of(state, op, tool);
    let subject = place == 0;
    let first = state
        .rows
        .iter()
        .find(|r| r.op == op)
        .is_some_and(|r| std::ptr::eq(r, row));

    let label = if first { op } else { "" };

    let colour = colour_at(place);
    let weight = if subject { BOLD } else { "" };

    let Some((wall, stddev)) = row.stats() else {
        return format!("  {DIM}{label:<11} {tool:<13}       —{RESET}");
    };

    // still running: a spinner rather than a bar, so a row
    // that is working looks like the status line does
    if row.walls.len() < state.repeat {
        let frame = SPINNER[(state.started.elapsed().as_secs_f64() / TICK.as_secs_f64()) as usize
            % SPINNER.len()];

        return format!(
            "  {label:<11} {weight}{colour}{tool:<13}{RESET} {colour}{frame}{RESET} \
             {DIM}{}/{}{RESET}",
            row.walls.len(),
            state.repeat,
        );
    }

    let ratio = subject_of(state, op).map(|s| wall / s);

    let ratio_text = match ratio {
        Some(r) if (r - 1.0).abs() < 0.005 => format!("{DIM}1.00x{RESET}"),
        Some(r) if r > 1.0 => format!("{RED}{r:.2}x{RESET}"),
        Some(r) => format!("{GREEN}{r:.2}x{RESET}"),
        None => "—".to_string(),
    };

    format!(
        "  {label:<11} {weight}{colour}{tool:<13} {:>11} {DIM}{:>9}{RESET}{weight}{colour} {:>10}{RESET} {:>8}",
        ms(wall),
        ms(stddev),
        // ru_maxrss is KiB on linux, so this is MiB and
        // says so: the difference reaches 5% by a gigabyte
        format!("{:.1} MiB", row.max_rss_kb as f64 / 1024.0),
        ratio_text,
    )
}

/// Rewrite the whole board in place.
//
// the cursor goes back over what was painted and every line
// is cleared before it is written, because a shorter line
// would otherwise leave the tail of a longer one behind
fn paint(state: &mut State) {
    let mut out = std::io::stdout().lock();

    if state.painted > 0 {
        let _ = write!(out, "\x1b[{}A", state.painted);
    }

    let mut lines = vec![
        format!("\x1b[2K {DIM}{}{RESET}", state.header),
        "\x1b[2K".to_string(),
        format!(
            "\x1b[2K  {BOLD}{:<11} {:<13} {:>11} {:>9} {:>10} {:>8}{RESET}",
            "operation",
            "tool",
            "wall",
            "±",
            "rss",
            format!("vs {}", state.subject)
        ),
        format!("\x1b[2K  {DIM}{}{RESET}", "─".repeat(67)),
    ];

    let ops: Vec<&str> = {
        let mut seen = Vec::new();
        for r in &state.rows {
            if !seen.contains(&r.op) {
                seen.push(r.op);
            }
        }
        seen
    };

    for (n, op) in ops.iter().enumerate() {
        // a rule between groups, so the eye reads one
        // operation's tools together rather than a column
        if n > 0 {
            lines.push(format!("\x1b[2K  {DIM}{}{RESET}", "┄".repeat(67)));
        }

        for tool in state
            .rows
            .iter()
            .filter(|r| r.op == *op)
            .map(|r| r.tool)
            .collect::<Vec<_>>()
        {
            lines.push(format!("\x1b[2K{}", line_for(state, op, tool)));
        }
    }

    lines.push("\x1b[2K".to_string());
    for line in status(state) {
        lines.push(format!("\x1b[2K{line}"));
    }

    for line in &lines {
        let _ = writeln!(out, "{line}");
    }
    let _ = out.flush();

    state.painted = lines.len();
}

/// The line under the board: what is running, how long it has been, and how
/// far through the whole run this is.
//
// the command itself rather than a step name, because what
// is worth watching is which tool is on the clock right now
fn status(state: &State) -> Vec<String> {
    let done = state.rows.iter().map(|r| r.walls.len()).sum::<usize>();
    let want = state.rows.len() * state.repeat;
    let total = state.started.elapsed().as_secs_f64();

    let counter = format!("  {DIM}{done}/{want} runs · {total:.1}s{RESET}");

    let Some(running) = &state.running else {
        return vec![counter];
    };

    let frame = SPINNER[(total / TICK.as_secs_f64()) as usize % SPINNER.len()];
    let held = running.since.elapsed().as_secs_f64();

    let shown = short(&running.line);

    vec![
        format!(
            "  {MAGENTA}{frame}{RESET} {BOLD}{:<14}{RESET}{DIM}{held:>5.1}s{RESET}  {CYAN}{shown}{RESET}",
            running.label,
        ),
        counter,
    ]
}

/// A command as something worth watching go past: no redirects, and every path
/// down to its last segment.
//
// the line michi records is meant to be pasteable, so it
// carries absolute paths and both redirects. that is the
// right thing for the table and the wrong thing for a line
// that changes ten times a second
fn short(line: &str) -> String {
    let body = line
        .split(" > ")
        .next()
        .unwrap_or(line)
        .split(" 2> ")
        .next()
        .unwrap_or(line);

    let words: Vec<String> = body
        .split_whitespace()
        .map(|word| match word.rsplit_once('/') {
            Some((_, last)) if !last.is_empty() => last.to_string(),
            _ => word.to_string(),
        })
        .collect();

    let shown = words.join(" ");
    let width = 76;

    if shown.chars().count() > width {
        format!("{}…", shown.chars().take(width).collect::<String>())
    } else {
        shown
    }
}

/// A line with its colour taken back out, for a run that is not going to a
/// terminal.
//
// stripped after formatting rather than never added:
// padding is applied to the values, so the columns line up
// either way and there is one code path to read
fn plain(line: &str) -> String {
    let mut out = String::with_capacity(line.len());
    let mut chars = line.chars();

    while let Some(c) = chars.next() {
        if c != '\x1b' {
            out.push(c);
            continue;
        }

        // an escape runs to the letter that ends it
        for c in chars.by_ref() {
            if c.is_ascii_alphabetic() {
                break;
            }
        }
    }

    out
}
