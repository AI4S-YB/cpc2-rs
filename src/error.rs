use std::error::Error;
use std::fmt::{self, Display, Formatter};
use std::io;

#[derive(Debug)]
pub enum Cpc2Error {
    Io(io::Error),
    Parse(String),
    Args(String),
}

impl Display for Cpc2Error {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        match self {
            Self::Io(err) => write!(f, "{err}"),
            Self::Parse(message) | Self::Args(message) => write!(f, "{message}"),
        }
    }
}

impl Error for Cpc2Error {
    fn source(&self) -> Option<&(dyn Error + 'static)> {
        match self {
            Self::Io(err) => Some(err),
            Self::Parse(_) | Self::Args(_) => None,
        }
    }
}

impl From<io::Error> for Cpc2Error {
    fn from(value: io::Error) -> Self {
        Self::Io(value)
    }
}
