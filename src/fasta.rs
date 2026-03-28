use std::io::{BufRead, Lines};

use crate::error::Cpc2Error;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FastaRecord {
    pub id: String,
    pub sequence: String,
}

impl FastaRecord {
    fn new(header: String, sequence: String) -> Result<Self, Cpc2Error> {
        let id = header
            .split_whitespace()
            .next()
            .ok_or_else(|| Cpc2Error::Parse("encountered an empty FASTA header".to_string()))?;

        Ok(Self {
            id: id.to_string(),
            sequence,
        })
    }
}

pub struct FastaReader<R: BufRead> {
    lines: Lines<R>,
    pending_header: Option<String>,
    finished: bool,
}

impl<R: BufRead> FastaReader<R> {
    pub fn new(reader: R) -> Self {
        Self {
            lines: reader.lines(),
            pending_header: None,
            finished: false,
        }
    }
}

impl<R: BufRead> Iterator for FastaReader<R> {
    type Item = Result<FastaRecord, Cpc2Error>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.finished {
            return None;
        }

        let header = match self.pending_header.take() {
            Some(header) => header,
            None => loop {
                let line = match self.lines.next()? {
                    Ok(line) => line,
                    Err(err) => return Some(Err(err.into())),
                };

                let trimmed = line.trim();
                if trimmed.is_empty() {
                    continue;
                }
                if let Some(header) = trimmed.strip_prefix('>') {
                    break header.trim().to_string();
                }
                return Some(Err(Cpc2Error::Parse(format!(
                    "FASTA record must start with '>', got: {trimmed}"
                ))));
            },
        };

        let mut sequence = String::new();

        while let Some(line) = self.lines.next() {
            let line = match line {
                Ok(line) => line,
                Err(err) => return Some(Err(err.into())),
            };
            let trimmed = line.trim();
            if trimmed.is_empty() {
                continue;
            }
            if let Some(next_header) = trimmed.strip_prefix('>') {
                self.pending_header = Some(next_header.trim().to_string());
                return Some(FastaRecord::new(header, sequence));
            }
            sequence.push_str(trimmed);
        }

        self.finished = true;
        Some(FastaRecord::new(header, sequence))
    }
}

#[cfg(test)]
mod tests {
    use std::io::Cursor;

    use super::FastaReader;

    #[test]
    fn parses_multiline_fasta() {
        let input = b">seq1 description\nATG\nTAA\n\n>seq2\nCCCGGG\n";
        let records: Vec<_> = FastaReader::new(Cursor::new(&input[..]))
            .map(|record| record.expect("valid record"))
            .collect();

        assert_eq!(records.len(), 2);
        assert_eq!(records[0].id, "seq1");
        assert_eq!(records[0].sequence, "ATGTAA");
        assert_eq!(records[1].id, "seq2");
        assert_eq!(records[1].sequence, "CCCGGG");
    }
}
