pub fn translate_dna(sequence: &str) -> String {
    let mut peptide = String::with_capacity(sequence.len() / 3);
    let mut index = 0;

    while index + 3 <= sequence.len() {
        peptide.push(codon_to_amino_acid(&sequence[index..index + 3]));
        index += 3;
    }

    peptide
}

pub fn sanitized_peptide(sequence: &str) -> String {
    sequence
        .chars()
        .filter(|aa| !matches!(aa, 'X' | 'B' | 'Z' | 'J' | 'U'))
        .collect()
}

pub fn trimmed_peptide(sequence: &str) -> String {
    sanitized_peptide(sequence).trim_matches('*').to_string()
}

pub fn isoelectric_point(sequence: &str) -> f64 {
    if sequence.is_empty() {
        return 0.0;
    }

    let mut low = 4.05;
    let mut high = 12.0;
    let mut ph = 7.775;

    while high - low > 0.0001 {
        let charge = charge_at_ph(sequence, ph);
        if charge > 0.0 {
            low = ph;
        } else {
            high = ph;
        }
        ph = (low + high) / 2.0;
    }

    ph
}

fn charge_at_ph(sequence: &str, ph: f64) -> f64 {
    let mut counts = ChargedAaCounts::default();
    for aa in sequence.chars() {
        match aa {
            'K' => counts.k += 1.0,
            'R' => counts.r += 1.0,
            'H' => counts.h += 1.0,
            'D' => counts.d += 1.0,
            'E' => counts.e += 1.0,
            'C' => counts.c += 1.0,
            'Y' => counts.y += 1.0,
            _ => {}
        }
    }

    let first = sequence.chars().next().unwrap_or('X');
    let last = sequence.chars().last().unwrap_or('X');

    let positive_charge = positive_partial(ph, n_terminal_pk(first))
        + counts.k * positive_partial(ph, 10.0)
        + counts.r * positive_partial(ph, 12.0)
        + counts.h * positive_partial(ph, 5.98);

    let negative_charge = negative_partial(ph, c_terminal_pk(last))
        + counts.d * negative_partial(ph, 4.05)
        + counts.e * negative_partial(ph, 4.45)
        + counts.c * negative_partial(ph, 9.0)
        + counts.y * negative_partial(ph, 10.0);

    positive_charge - negative_charge
}

fn positive_partial(ph: f64, pk: f64) -> f64 {
    1.0 / (10f64.powf(ph - pk) + 1.0)
}

fn negative_partial(ph: f64, pk: f64) -> f64 {
    1.0 / (10f64.powf(pk - ph) + 1.0)
}

fn n_terminal_pk(aa: char) -> f64 {
    match aa {
        'A' => 7.59,
        'M' => 7.0,
        'S' => 6.93,
        'P' => 8.36,
        'T' => 6.82,
        'V' => 7.44,
        'E' => 7.7,
        _ => 7.5,
    }
}

fn c_terminal_pk(aa: char) -> f64 {
    match aa {
        'D' => 4.55,
        'E' => 4.75,
        _ => 3.55,
    }
}

#[derive(Default)]
struct ChargedAaCounts {
    k: f64,
    r: f64,
    h: f64,
    d: f64,
    e: f64,
    c: f64,
    y: f64,
}

fn codon_to_amino_acid(codon: &str) -> char {
    match codon {
        "TTT" | "TTC" => 'F',
        "TTA" | "TTG" | "CTT" | "CTC" | "CTA" | "CTG" => 'L',
        "ATT" | "ATC" | "ATA" => 'I',
        "ATG" => 'M',
        "GTT" | "GTC" | "GTA" | "GTG" => 'V',
        "TCT" | "TCC" | "TCA" | "TCG" | "AGT" | "AGC" => 'S',
        "CCT" | "CCC" | "CCA" | "CCG" => 'P',
        "ACT" | "ACC" | "ACA" | "ACG" => 'T',
        "GCT" | "GCC" | "GCA" | "GCG" => 'A',
        "TAT" | "TAC" => 'Y',
        "TAA" | "TAG" | "TGA" => '*',
        "CAT" | "CAC" => 'H',
        "CAA" | "CAG" => 'Q',
        "AAT" | "AAC" => 'N',
        "AAA" | "AAG" => 'K',
        "GAT" | "GAC" => 'D',
        "GAA" | "GAG" => 'E',
        "TGT" | "TGC" => 'C',
        "TGG" => 'W',
        "CGT" | "CGC" | "CGA" | "CGG" | "AGA" | "AGG" => 'R',
        "GGT" | "GGC" | "GGA" | "GGG" => 'G',
        _ => 'X',
    }
}

#[cfg(test)]
mod tests {
    use super::{isoelectric_point, translate_dna, trimmed_peptide};

    #[test]
    fn translates_standard_codons() {
        assert_eq!(translate_dna("ATGGGCTAA"), "MG*");
    }

    #[test]
    fn trims_terminal_stop() {
        assert_eq!(trimmed_peptide("MG*"), "MG");
    }

    #[test]
    fn matches_biopython_examples() {
        let high = isoelectric_point("INGAR");
        let low = isoelectric_point("PETER");
        assert!((high - 9.75).abs() < 0.01);
        assert!((low - 4.53).abs() < 0.01);
    }
}
