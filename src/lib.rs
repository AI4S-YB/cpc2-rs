pub mod error;
pub mod fasta;
pub mod fickett;
pub mod orf;
pub mod protein;
pub mod svm;

use std::fmt::{self, Display, Formatter};
use std::io::{BufRead, Write};

use error::Cpc2Error;
use fasta::FastaReader;
use fickett::Fickett;
use orf::longest_orf;
use protein::{isoelectric_point, translate_dna, trimmed_peptide};
use svm::embedded_predictor;

#[derive(Debug, Clone, Copy)]
pub struct AnalyzeOptions {
    pub check_reverse: bool,
    pub include_orf_start: bool,
    pub include_peptide: bool,
}

#[derive(Debug, Clone)]
pub struct TranscriptPrediction {
    pub id: String,
    pub transcript_length: usize,
    pub peptide_length: usize,
    pub fickett_score: f64,
    pub isoelectric_point: f64,
    pub orf_integrity: i8,
    pub orf_start: usize,
    pub putative_peptide: String,
    pub coding_probability: f64,
    pub label: CodingLabel,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CodingLabel {
    Coding,
    Noncoding,
}

impl Display for CodingLabel {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        match self {
            Self::Coding => write!(f, "coding"),
            Self::Noncoding => write!(f, "noncoding"),
        }
    }
}

pub fn analyze_reader<R: BufRead>(
    reader: R,
    options: &AnalyzeOptions,
) -> Result<Vec<TranscriptPrediction>, Cpc2Error> {
    let predictor = embedded_predictor();
    let mut predictions = Vec::new();

    for record in FastaReader::new(reader) {
        let record = record?;
        let normalized_sequence = normalize_sequence(&record.sequence);
        let transcript_length = normalized_sequence.len();
        let longest = longest_orf(&normalized_sequence, options.check_reverse);
        let translated = translate_dna(&longest.sequence);
        let peptide_length = translated.len();
        let cleaned_peptide = trimmed_peptide(&translated);
        let fickett_score = Fickett::score(&normalized_sequence);

        let (orf_start, orf_integrity, isoelectric_point) =
            if peptide_length > 0 && !cleaned_peptide.is_empty() {
                (
                    longest.start + 1,
                    longest.integrity,
                    isoelectric_point(&cleaned_peptide),
                )
            } else {
                (0, -1, 0.0)
            };

        let prediction = predictor.predict([
            peptide_length as f64,
            fickett_score,
            isoelectric_point,
            orf_integrity as f64,
        ]);

        predictions.push(TranscriptPrediction {
            id: record.id,
            transcript_length,
            peptide_length,
            fickett_score,
            isoelectric_point,
            orf_integrity,
            orf_start,
            putative_peptide: if peptide_length > 0 {
                cleaned_peptide
            } else {
                "non".to_string()
            },
            coding_probability: prediction.coding_probability,
            label: if prediction.predicted_label == 1 {
                CodingLabel::Coding
            } else {
                CodingLabel::Noncoding
            },
        });
    }

    Ok(predictions)
}

pub fn write_results<W: Write>(
    writer: &mut W,
    results: &[TranscriptPrediction],
    options: &AnalyzeOptions,
) -> Result<(), Cpc2Error> {
    let mut header = vec![
        "#ID".to_string(),
        "transcript_length".to_string(),
        "peptide_length".to_string(),
        "Fickett_score".to_string(),
        "pI".to_string(),
        "ORF_integrity".to_string(),
    ];

    if options.include_orf_start {
        header.push("ORF_Start".to_string());
    }
    if options.include_peptide {
        header.push("putative_peptide".to_string());
    }
    header.push("coding_probability".to_string());
    header.push("label".to_string());

    writeln!(writer, "{}", header.join("\t"))?;

    for result in results {
        let mut columns = vec![
            result.id.clone(),
            result.transcript_length.to_string(),
            result.peptide_length.to_string(),
            format_float(result.fickett_score),
            format_float(result.isoelectric_point),
            result.orf_integrity.to_string(),
        ];

        if options.include_orf_start {
            columns.push(result.orf_start.to_string());
        }
        if options.include_peptide {
            columns.push(result.putative_peptide.clone());
        }
        columns.push(format_float(result.coding_probability));
        columns.push(result.label.to_string());

        writeln!(writer, "{}", columns.join("\t"))?;
    }

    Ok(())
}

fn normalize_sequence(sequence: &str) -> String {
    sequence
        .trim()
        .chars()
        .map(|base| match base.to_ascii_uppercase() {
            'U' => 'T',
            other => other,
        })
        .collect()
}

fn format_float(value: f64) -> String {
    let mut output = format!("{value:.10}");
    while output.contains('.') && output.ends_with('0') {
        output.pop();
    }
    if output.ends_with('.') {
        output.pop();
    }
    if !output.contains('.') {
        output.push_str(".0");
    }
    output
}

#[cfg(test)]
mod tests {
    use std::io::Cursor;

    use super::{analyze_reader, AnalyzeOptions, CodingLabel};

    #[test]
    fn analyzes_small_example() {
        let input = b">seq1\nATGAAATAG\n";
        let results = analyze_reader(
            Cursor::new(&input[..]),
            &AnalyzeOptions {
                check_reverse: false,
                include_orf_start: true,
                include_peptide: true,
            },
        )
        .expect("analysis should succeed");

        assert_eq!(results.len(), 1);
        assert_eq!(results[0].id, "seq1");
        assert_eq!(results[0].peptide_length, 3);
        assert_eq!(results[0].orf_start, 1);
        assert_eq!(results[0].putative_peptide, "MK");
        assert!(matches!(
            results[0].label,
            CodingLabel::Coding | CodingLabel::Noncoding
        ));
    }

    #[test]
    fn matches_example_forward_mode() {
        let results = analyze_reader(
            Cursor::new(include_str!("../data/example.fa").as_bytes()),
            &AnalyzeOptions {
                check_reverse: false,
                include_orf_start: true,
                include_peptide: true,
            },
        )
        .expect("analysis should succeed");

        assert_eq!(results.len(), 2);

        assert_eq!(results[0].id, "AF282387");
        assert_eq!(results[0].peptide_length, 176);
        assert_eq!(results[0].orf_start, 1);
        assert_eq!(results[0].label, CodingLabel::Coding);
        assert!((results[0].coding_probability - 0.997542).abs() < 1e-4);

        assert_eq!(results[1].id, "Tsix_mus");
        assert_eq!(results[1].peptide_length, 70);
        assert_eq!(results[1].orf_start, 105);
        assert_eq!(
            results[1].putative_peptide,
            "MAAVTGRTVAACYNQKTTSGFSLCPKKVVKCANHTLSAGSRASRALGPAPHSLKPSVQRLCQAQSRKIR"
        );
        assert_eq!(results[1].label, CodingLabel::Noncoding);
        assert!((results[1].coding_probability - 0.0447345).abs() < 1e-4);
    }

    #[test]
    fn matches_example_reverse_mode() {
        let results = analyze_reader(
            Cursor::new(include_str!("../data/example.fa").as_bytes()),
            &AnalyzeOptions {
                check_reverse: true,
                include_orf_start: true,
                include_peptide: true,
            },
        )
        .expect("analysis should succeed");

        assert_eq!(results.len(), 2);
        assert_eq!(results[1].id, "Tsix_mus");
        assert_eq!(results[1].peptide_length, 80);
        assert_eq!(results[1].orf_start, 1643);
        assert_eq!(
            results[1].putative_peptide,
            "MSELQSLWPLLFWSLRLQRRGSVKGVQRLAPAMFARFPWMCGSSVVSLHLRSFGGTFLVPLPPSLMAYLRKHIKIPRDF"
        );
        assert_eq!(results[1].label, CodingLabel::Noncoding);
        assert!((results[1].coding_probability - 0.0521955).abs() < 1e-4);
    }
}
