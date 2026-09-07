//! Reading a file at an absolute offset, on the platforms that can.
//!
//! `std` spells this two ways and neither is portable, so the crate spells it
//! once here rather than once per format that reads at an offset. The Windows
//! arm has never been run: it typechecks under `cargo clippy --target
//! x86_64-pc-windows-msvc` and that is the whole of the evidence.
//!
//! Both arms take `&File` rather than `&mut File`, which is the whole point:
//! [`Indexable::get`](crate::collection::Indexable::get) takes `&self`, so a
//! disk-backed collection reading through `Seek` would need a `RefCell` or a
//! `Mutex` to reach its handle. A positioned read needs neither, and lets
//! several threads read one collection at once.

use std::fs::File;
use std::io;

/// Fill `buf` from `offset`, leaving the file's cursor alone.
#[cfg(unix)]
pub fn read_exact_at(file: &File, buf: &mut [u8], offset: u64) -> io::Result<()> {
    use std::os::unix::fs::FileExt;

    file.read_exact_at(buf, offset)
}

/// Fill `buf` from `offset`, leaving the file's cursor alone.
#[cfg(windows)]
pub fn read_exact_at(file: &File, mut buf: &mut [u8], mut offset: u64) -> io::Result<()> {
    use std::os::windows::fs::FileExt;

    // seek_read may return short and has no read_exact
    // counterpart, so the filling loop is written out
    while !buf.is_empty() {
        // **note: unlike pread this moves the file cursor.
        //         nothing in the crate reads that cursor
        //         after open, which is the only reason the
        //         two arms are interchangeable
        match file.seek_read(buf, offset) {
            Ok(0) => break,
            Ok(n) => {
                buf = &mut buf[n..];
                offset += n as u64;
            }
            Err(e) if e.kind() == io::ErrorKind::Interrupted => continue,
            Err(e) => return Err(e),
        }
    }

    if !buf.is_empty() {
        // the wording std's read_exact uses, since callers
        // report it and one test matches on it
        return Err(io::Error::new(
            io::ErrorKind::UnexpectedEof,
            "failed to fill whole buffer",
        ));
    }

    Ok(())
}
