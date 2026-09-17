//! Drop a file's pages from the page cache, so the next read of it reaches the
//! disk.
//!
//! Every benchmark here reads the same gigabyte over and over, and a machine
//! with room to keep it resident will do exactly that. Left alone the suite
//! measures memcpy out of the page cache, which flatters anything whose cost is
//! reading and hides everything a locality difference would show.

use std::fs::File;
use std::os::fd::AsRawFd;

fn main() {
    for path in std::env::args().skip(1) {
        let Ok(file) = File::open(&path) else {
            // a run that resets before every repeat asks for
            // files that are not there yet, and a missing
            // one is already evicted
            continue;
        };

        // POSIX_FADV_DONTNEED drops clean pages and leaves
        // dirty ones, so a fixture written moments ago would
        // survive its own eviction without this
        let _ = file.sync_all();

        // SAFETY: the fd is open for the length of the call,
        //         and posix_fadvise reads no memory through
        //         it -- 0, 0 names the whole file
        let rc = unsafe { libc::posix_fadvise(file.as_raw_fd(), 0, 0, libc::POSIX_FADV_DONTNEED) };

        assert!(rc == 0, "posix_fadvise({path}) failed with {rc}");
    }
}
