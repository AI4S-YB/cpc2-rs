#[derive(Debug, Clone, PartialEq, Eq)]
pub struct OrfResult {
    pub sequence: String,
    pub start: usize,
    pub strand: char,
    pub integrity: i8,
}

#[derive(Debug, Clone, Copy)]
struct BestOrf {
    direction: char,
    start: usize,
    end: usize,
    length: usize,
    integrity: i8,
}

impl Default for BestOrf {
    fn default() -> Self {
        Self {
            direction: '+',
            start: 0,
            end: 0,
            length: 0,
            integrity: 0,
        }
    }
}

pub fn longest_orf(sequence: &str, check_reverse: bool) -> OrfResult {
    let mut best = BestOrf::default();
    scan_frames(sequence.as_bytes(), '+', &mut best);

    let reverse = if check_reverse {
        let reverse = reverse_complement(sequence);
        scan_frames(reverse.as_bytes(), '-', &mut best);
        Some(reverse)
    } else {
        None
    };

    let source = if best.direction == '-' {
        reverse.as_deref().unwrap_or(sequence)
    } else {
        sequence
    };

    OrfResult {
        sequence: source[best.start..best.end].to_string(),
        start: best.start,
        strand: best.direction,
        integrity: best.integrity,
    }
}

fn scan_frames(sequence: &[u8], direction: char, best: &mut BestOrf) {
    for frame in 0..3 {
        find_longest_in_frame(sequence, frame, direction, best);
    }
}

fn find_longest_in_frame(sequence: &[u8], frame: usize, direction: char, best: &mut BestOrf) {
    let len = sequence.len();
    let mut coordinate = frame;

    while coordinate + 3 <= len {
        let codon = &sequence[coordinate..coordinate + 3];
        if is_start_codon(codon) && !is_stop_codon(codon) {
            let orf_start = coordinate;
            let mut integrity = -1;
            let mut last_index = coordinate;

            loop {
                coordinate += 3;
                if coordinate + 3 > len {
                    break;
                }

                last_index = coordinate;
                let codon = &sequence[coordinate..coordinate + 3];
                if is_stop_codon(codon) {
                    integrity = 1;
                    break;
                }
            }

            let orf_end = if integrity == 1 {
                coordinate + 3
            } else {
                last_index + 3
            };
            let length = orf_end - orf_start;

            if length > best.length || (length == best.length && orf_start < best.start) {
                best.direction = direction;
                best.start = orf_start;
                best.end = orf_end;
                best.length = length;
                best.integrity = integrity;
            }

            if integrity == 1 {
                coordinate += 3;
            }
            continue;
        }

        coordinate += 3;
    }
}

fn is_start_codon(codon: &[u8]) -> bool {
    codon == b"ATG"
}

fn is_stop_codon(codon: &[u8]) -> bool {
    codon == b"TAA" || codon == b"TAG" || codon == b"TGA"
}

fn reverse_complement(sequence: &str) -> String {
    sequence
        .chars()
        .rev()
        .map(|base| match base.to_ascii_uppercase() {
            'A' => 'T',
            'T' => 'A',
            'U' => 'A',
            'C' => 'G',
            'G' => 'C',
            'N' => 'N',
            'X' => 'X',
            _ => 'N',
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::longest_orf;

    #[test]
    fn finds_complete_orf() {
        let orf = longest_orf("CCCATGAAATAGGG", false);
        assert_eq!(orf.sequence, "ATGAAATAG");
        assert_eq!(orf.start, 3);
        assert_eq!(orf.integrity, 1);
        assert_eq!(orf.strand, '+');
    }
}
