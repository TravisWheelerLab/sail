//! The crate's error type.
//!
//! Two variants, for the two things that go wrong here: a record fails its
//! own invariants, or writing one out fails at the sink.

use std::fmt;

/// A [`Result`](std::result::Result) with this crate's error type.
pub type Result<T> = std::result::Result<T, Error>;

/// What went wrong. `#[non_exhaustive]`, so match with a `_` arm.
#[derive(Debug)]
#[non_exhaustive]
pub enum ErrorKind {
    /// A record does not satisfy its format's rules.
    Parse { detail: String },
    /// Writing a record out failed.
    Io { source: std::io::Error },
}

/// An error from building or writing a record.
#[derive(Debug)]
pub struct Error {
    // held inline rather than boxed: no Result in this
    // crate sits on a per-line path, and a two-variant enum
    // is cheaper to carry than a pointer chase to reach it
    kind: ErrorKind,
}

impl Error {
    pub fn parse(detail: impl Into<String>) -> Self {
        Error {
            kind: ErrorKind::Parse {
                detail: detail.into(),
            },
        }
    }

    /// What kind of failure this is, coarsely enough to branch on.
    pub fn kind(&self) -> &ErrorKind {
        &self.kind
    }
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match &self.kind {
            ErrorKind::Parse { detail } => write!(f, "{detail}"),
            ErrorKind::Io { source } => write!(f, "{source}"),
        }
    }
}

impl std::error::Error for Error {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match &self.kind {
            ErrorKind::Io { source } => Some(source),
            ErrorKind::Parse { .. } => None,
        }
    }
}

impl From<std::io::Error> for Error {
    fn from(source: std::io::Error) -> Self {
        Error {
            kind: ErrorKind::Io { source },
        }
    }
}

/// Return an [`ErrorKind::Parse`] error, formatted.
#[macro_export]
macro_rules! bail_parse {
    ($($arg:tt)*) => {
        return Err($crate::Error::parse(format!($($arg)*)))
    };
}

pub use bail_parse;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn an_io_failure_keeps_its_cause_reachable() {
        use std::error::Error as _;

        let err: Error = std::io::Error::new(std::io::ErrorKind::BrokenPipe, "nope").into();
        assert!(matches!(err.kind(), ErrorKind::Io { .. }));
        assert_eq!(err.to_string(), "nope");
        assert!(err.source().is_some(), "the io error stays reachable");
    }

    #[test]
    fn a_parse_failure_says_only_what_is_wrong() {
        let err = Error::parse("rows are 4 and 6 columns wide");
        assert!(matches!(err.kind(), ErrorKind::Parse { .. }));
        assert_eq!(err.to_string(), "rows are 4 and 6 columns wide");
        assert!(std::error::Error::source(&err).is_none());
    }
}
