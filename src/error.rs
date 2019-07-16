use std::fmt;

#[derive(Debug)]
pub struct Error {
    message: String,
}

impl Error {
    pub fn new<S>(message: S) -> Error
    where
        String: From<S>,
    {
        Error {
            message: message.into(),
        }
    }
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}", self.message)
    }
}

impl std::error::Error for Error {}

#[macro_export]
macro_rules! fail {
    ($e:expr) => {
        return Err(Error::new($e).into());
    };
    ($fmt:expr, $($arg:tt)+) => {
        return Err(Error::new(format!($fmt, $($arg)+)).into());
    };
}
