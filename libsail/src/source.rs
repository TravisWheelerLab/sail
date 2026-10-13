//! Where an indexed collection's bytes come from.
//!
//! [`FileSource`] reads each range when asked; `MmapSource`, with the `mmap`
//! feature, has the file mapped and hands back a slice of it. Which one a
//! collection holds is a type parameter, so a caller picks the trade and pays
//! for nothing else.

use std::borrow::Cow;
use std::fs::File;
use std::io;
use std::path::Path;

/// Bytes addressable at an offset, without a cursor.
pub trait Source {
    /// The `len` bytes at `start`.
    //
    // **note: `&self`, not `&mut self`. Indexable::get takes
    //         &self, so a Read or Seek bound would put a
    //         RefCell or a Mutex behind every disk-backed
    //         collection -- the whole reason pread.rs
    //         exists. see src/pread.rs
    //
    // Cow, not a caller-supplied buffer: a mapped source
    // borrows and copies nothing, which is most of what it
    // is for. filling &mut [u8] would throw that away
    fn range(&self, start: u64, len: usize) -> io::Result<Cow<'_, [u8]>>;
}

// ---

/// A file each range is read from when asked. Costs a syscall and an allocation
/// per range, and keeps only the handle resident.
pub struct FileSource {
    file: File,
}

impl FileSource {
    pub fn open(path: impl AsRef<Path>) -> io::Result<Self> {
        Ok(FileSource {
            file: File::open(path)?,
        })
    }

    /// The handle, so a collection can index before asking for a range.
    pub fn file(&self) -> &File {
        &self.file
    }
}

impl Source for FileSource {
    fn range(&self, start: u64, len: usize) -> io::Result<Cow<'_, [u8]>> {
        let mut buf = vec![0u8; len];
        crate::pread::read_exact_at(&self.file, &mut buf, start)?;

        Ok(Cow::Owned(buf))
    }
}

// ---

/// A file mapped into memory, so a range is a slice of it. About 1.66x faster
/// per range than [`FileSource`], at the cost of resident pages.
///
/// # Safety
///
/// Constructing one is unsafe; see [`open`](Self::open).
#[cfg(feature = "mmap")]
pub struct MmapSource {
    // the resident cost is page granularity and not bytes: a
    // 400 byte record costs a whole page, so 100,000 records
    // drawn at random over a 1 GiB file leave ~78% of it
    // resident. those pages are clean and reclaimable, and a
    // FileSource puts the same file in the page cache
    // without being billed for it -- the process is billed
    // here. measured by benchmarks/bin/mmap.rs
    map: memmap2::Mmap,
}

#[cfg(feature = "mmap")]
impl MmapSource {
    /// Map the file at `path`.
    ///
    /// # Safety
    ///
    /// The file must not be modified or truncated by anything, this process
    /// included, while the source is alive. A truncation under an open map
    /// raises `SIGBUS`, which aborts the process rather than returning an
    /// [`Error`](crate::Error) a caller could handle.
    pub unsafe fn open(path: impl AsRef<Path>) -> io::Result<Self> {
        let file = File::open(path)?;

        // SAFETY: nothing may modify or truncate the file
        //         while the map is alive. that is this
        //         function's own unsafe contract, which is
        //         the only place a caller can honour it
        let map = unsafe { memmap2::Mmap::map(&file)? };

        Ok(MmapSource { map })
    }

    /// The mapped bytes, so a collection can index before asking for a range.
    pub fn bytes(&self) -> &[u8] {
        &self.map
    }
}

#[cfg(feature = "mmap")]
impl Source for MmapSource {
    fn range(&self, start: u64, len: usize) -> io::Result<Cow<'_, [u8]>> {
        // a checked slice, so a map shorter than its offsets
        // fails the way a short file does rather than
        // panicking. only a truncation *after* mapping is
        // the SIGBUS that open's contract covers
        let start = usize::try_from(start).map_err(|_| short())?;
        let end = start.checked_add(len).ok_or_else(short)?;

        self.map
            .get(start..end)
            .map(Cow::Borrowed)
            .ok_or_else(short)
    }
}

/// The error a range past the end reports.
#[cfg(feature = "mmap")]
fn short() -> io::Error {
    // std's own wording for a read_exact that ran out, which
    // pread.rs's Windows arm reproduces too: tests match on
    // the string, and a source failing differently would
    // read as a different bug
    io::Error::new(io::ErrorKind::UnexpectedEof, "failed to fill whole buffer")
}
