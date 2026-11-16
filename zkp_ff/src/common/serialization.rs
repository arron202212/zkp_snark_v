use std::io::{self, Read, Write};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SerializationMode {
    Text,
    Binary,
}

/// 单个元素的序列化接口。
///
/// write_payload / read_payload 只负责元素自身的内容，
/// 容器负责元素之间的分隔符。
pub trait Serialize: Sized {
    fn write_payload<W: Write>(
        &self,
        writer: &mut W,
        mode: SerializationMode,
    ) -> io::Result<()>;

    fn read_payload<R: Read>(
        reader: &mut R,
        mode: SerializationMode,
    ) -> io::Result<Self>;
}


fn read_token<R: Read>(reader: &mut R) -> io::Result<String> {
    let mut token = Vec::new();
    let mut byte = [0u8; 1];

    loop {
        match reader.read(&mut byte)? {
            0 => break,
            _ if byte[0].is_ascii_whitespace() => {
                if token.is_empty() {
                    continue;
                }
                break;
            }
            _ => token.push(byte[0]),
        }
    }

    if token.is_empty() {
        return Err(io::Error::new(
            io::ErrorKind::UnexpectedEof,
            "unexpected end of input",
        ));
    }

    String::from_utf8(token).map_err(|e| {
        io::Error::new(io::ErrorKind::InvalidData, e)
    })
}

fn write_text<T: std::fmt::Display, W: Write>(
    writer: &mut W,
    value: T,
) -> io::Result<()> {
    write!(writer, "{}", value)
}

fn parse_text<T: std::str::FromStr, R: Read>(
    reader: &mut R,
) -> io::Result<T>
where
    T::Err: std::fmt::Display,
{
    let token = read_token(reader)?;

    token.parse::<T>().map_err(|e| {
        io::Error::new(io::ErrorKind::InvalidData, e.to_string())
    })
}


impl Serialize for u64 {
    fn write_payload<W: Write>(
        &self,
        writer: &mut W,
        mode: SerializationMode,
    ) -> io::Result<()> {
        match mode {
            SerializationMode::Text => {
                write_text(writer, *self)
            }
            SerializationMode::Binary => {
                writer.write_all(&self.to_le_bytes())
            }
        }
    }

    fn read_payload<R: Read>(
        reader: &mut R,
        mode: SerializationMode,
    ) -> io::Result<Self> {
        match mode {
            SerializationMode::Text => parse_text(reader),
            SerializationMode::Binary => {
                let mut bytes = [0u8; 8];
                reader.read_exact(&mut bytes)?;
                Ok(u64::from_le_bytes(bytes))
            }
        }
    }
}

pub fn output_bool<W: Write>(
    writer: &mut W,
    value: bool,
) -> io::Result<()> {
    writeln!(writer, "{}", if value { 1 } else { 0 })
}

pub fn input_bool<R: Read>(
    reader: &mut R,
) -> io::Result<bool> {
    let value: u64 = parse_text(reader)?;

    match value {
        0 => Ok(false),
        1 => Ok(true),
        _ => Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "boolean value must be 0 or 1",
        )),
    }
}

pub fn output_bool_vector<W: Write>(
    writer: &mut W,
    values: &[bool],
) -> io::Result<()> {
    writeln!(writer, "{}", values.len())?;

    for &value in values {
        output_bool(writer, value)?;
    }

    Ok(())
}

pub fn input_bool_vector<R: Read>(
    reader: &mut R,
) -> io::Result<Vec<bool>> {
    let len: usize = parse_text(reader)?;

    // 防止异常输入导致不合理的内存分配。
    const MAX_ELEMENTS: usize = 10_000_000;

    if len > MAX_ELEMENTS {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "boolean vector is too large",
        ));
    }

    let mut values = Vec::with_capacity(len);

    for _ in 0..len {
        values.push(input_bool(reader)?);
    }

    Ok(values)
}