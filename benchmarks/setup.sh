#!/usr/bin/env bash
#
# Fetch the tools `compare` measures sail against, into benchmarks/tools.
#
#   ./setup.sh           fetch whatever is missing
#   ./setup.sh --force    fetch everything again
#
# Nothing is installed outside this directory: no package manager, no
# ~/.local/bin, no PATH edits. src/cli.rs looks for these binaries in
# tools/bin before it looks at PATH, so a run picks up the pinned versions
# here even on a machine that has its own seqkit.
#
# The versions are pinned because they are part of the measurement. A number
# is only comparable against another number taken with the same peer, so
# bumping either version below is a change to the benchmark and not a
# maintenance chore.

set -euo pipefail

SEQKIT_VERSION=v2.13.0
EASEL_VERSION=easel-0.49

here=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
tools=$here/tools
bin=$tools/bin
src=$tools/src

force=false
if [[ ${1-} == --force ]]; then
    force=true
elif [[ $# -gt 0 ]]; then
    echo "usage: setup.sh [--force]" >&2
    exit 2
fi

say() { printf '==> %s\n' "$*"; }

# Fail now, naming the whole list, rather than partway through a build.
check_deps() {
    local missing=()
    for dep in "$@"; do
        command -v "$dep" >/dev/null || missing+=("$dep")
    done

    if (( ${#missing[@]} )); then
        echo "setup.sh needs ${missing[*]}, which $( (( ${#missing[@]} == 1 )) && echo is || echo are) not installed" >&2
        exit 1
    fi
}

# --- seqkit
#
# A release binary rather than a build from source, because seqkit is Go and
# this machine has no toolchain for it. The published .md5.txt is checked,
# which catches a truncated download; it is served by the same host as the
# tarball, so it is not evidence against a compromised release.

seqkit() {
    if ! $force && [[ -x $bin/seqkit ]] \
        && "$bin/seqkit" version 2>/dev/null | grep -qF "${SEQKIT_VERSION}"; then
        say "seqkit ${SEQKIT_VERSION} already in tools/bin"
        return
    fi

    check_deps curl tar md5sum

    local url=https://github.com/shenwei356/seqkit/releases/download/${SEQKIT_VERSION}/seqkit_linux_amd64.tar.gz
    local tarball=$src/seqkit_linux_amd64.tar.gz

    say "fetching seqkit ${SEQKIT_VERSION}"
    mkdir -p "$src" "$bin"
    curl -fsSL -o "$tarball" "$url"

    # the .md5.txt names the file, so compare only the digest
    local want got
    want=$(curl -fsSL "${url}.md5.txt" | awk '{print $1}')
    got=$(md5sum "$tarball" | awk '{print $1}')
    if [[ $want != "$got" ]]; then
        echo "seqkit tarball md5 is $got, the release says $want" >&2
        exit 1
    fi

    tar xzf "$tarball" -C "$bin" seqkit
    rm -f "$tarball"

    say "seqkit -> $("$bin"/seqkit version)"
}

# --- easel
#
# Built from source: the Eddy lab ships no binaries, and the miniapps are
# normally reached through a HMMER install. Cloning easel alone keeps this to
# one small repo.

easel() {
    if ! $force && [[ -x $bin/esl-seqstat ]] \
        && "$bin/esl-seqstat" -h 2>/dev/null | grep -qF "${EASEL_VERSION#easel-}"; then
        say "easel ${EASEL_VERSION} already in tools/bin"
        return
    fi

    check_deps git autoconf gcc make

    local repo=$src/easel

    say "fetching easel ${EASEL_VERSION}"
    mkdir -p "$src" "$bin"
    rm -rf "$repo"
    git clone --quiet --depth 1 --branch "$EASEL_VERSION" \
        https://github.com/EddyRivasLab/easel.git "$repo"

    # easel builds with a wall of upstream -Wdiscarded-qualifiers warnings, so
    # the build goes to a log and only a failure prints it. the log stays
    # either way, since a warning that turns out to matter is unreachable once
    # it scrolls past
    local log=$src/easel-build.log
    say "building easel (log at ${log#"$here"/})"
    if ! (
        cd "$repo"
        # the repo ships configure.ac and not configure
        autoconf
        ./configure --quiet
        make -j"$(nproc)"
    ) >"$log" 2>&1; then
        echo "easel build failed; last 40 lines of $log:" >&2
        tail -40 "$log" >&2
        exit 1
    fi

    # every miniapp, not just the three compare uses today: they are already
    # built by this point, and a new comparison should not need a setup change.
    # the name filter drops the .pl, .man and .in files sitting beside them
    find "$repo/miniapps" -maxdepth 1 -type f -executable -name 'esl-*' ! -name '*.*' \
        -exec cp {} "$bin/" \;

    say "easel -> $("$bin"/esl-seqstat -h | sed -n 's/^# \(Easel .*\)/\1/p')"
}

# ---

mkdir -p "$bin" "$src"
seqkit
easel

say "tools in $bin:"
ls "$bin" | sed 's/^/    /'
