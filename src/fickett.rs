pub struct Fickett;

impl Fickett {
    const POSITION_PARAMETER: [f64; 10] = [1.9, 1.8, 1.7, 1.6, 1.5, 1.4, 1.3, 1.2, 1.1, 0.0];
    const CONTENT_PARAMETER: [f64; 10] =
        [0.33, 0.31, 0.29, 0.27, 0.25, 0.23, 0.21, 0.19, 0.17, 0.0];

    pub fn score(dna: &str) -> f64 {
        if dna.len() < 2 {
            return 0.0;
        }

        let total_base = dna.len() as f64;
        let phase_0 = dna.as_bytes().iter().step_by(3);
        let phase_1 = dna.as_bytes().iter().skip(1).step_by(3);
        let phase_2 = dna.as_bytes().iter().skip(2).step_by(3);

        let phase_counts = [
            count_bases(phase_0),
            count_bases(phase_1),
            count_bases(phase_2),
        ];

        let a_content =
            (phase_counts[0].0 + phase_counts[1].0 + phase_counts[2].0) as f64 / total_base;
        let c_content =
            (phase_counts[0].1 + phase_counts[1].1 + phase_counts[2].1) as f64 / total_base;
        let g_content =
            (phase_counts[0].2 + phase_counts[1].2 + phase_counts[2].2) as f64 / total_base;
        let t_content =
            (phase_counts[0].3 + phase_counts[1].3 + phase_counts[2].3) as f64 / total_base;

        let a_position = position_ratio([phase_counts[0].0, phase_counts[1].0, phase_counts[2].0]);
        let c_position = position_ratio([phase_counts[0].1, phase_counts[1].1, phase_counts[2].1]);
        let g_position = position_ratio([phase_counts[0].2, phase_counts[1].2, phase_counts[2].2]);
        let t_position = position_ratio([phase_counts[0].3, phase_counts[1].3, phase_counts[2].3]);

        Self::content_probability(a_content, 'A')
            + Self::content_probability(c_content, 'C')
            + Self::content_probability(g_content, 'G')
            + Self::content_probability(t_content, 'T')
            + Self::position_probability(a_position, 'A')
            + Self::position_probability(c_position, 'C')
            + Self::position_probability(g_position, 'G')
            + Self::position_probability(t_position, 'T')
    }

    fn position_probability(value: f64, base: char) -> f64 {
        if value < 0.0 {
            return 0.0;
        }

        for (idx, threshold) in Self::POSITION_PARAMETER.iter().enumerate() {
            if value >= *threshold {
                return position_probabilities(base)[idx] * position_weight(base);
            }
        }

        0.0
    }

    fn content_probability(value: f64, base: char) -> f64 {
        if value < 0.0 {
            return 0.0;
        }

        for (idx, threshold) in Self::CONTENT_PARAMETER.iter().enumerate() {
            if value >= *threshold {
                return content_probabilities(base)[idx] * content_weight(base);
            }
        }

        0.0
    }
}

fn count_bases<'a>(iter: impl Iterator<Item = &'a u8>) -> (usize, usize, usize, usize) {
    let mut counts = (0, 0, 0, 0);
    for base in iter {
        match *base {
            b'A' => counts.0 += 1,
            b'C' => counts.1 += 1,
            b'G' => counts.2 += 1,
            b'T' => counts.3 += 1,
            _ => {}
        }
    }
    counts
}

fn position_ratio(values: [usize; 3]) -> f64 {
    let min = *values.iter().min().unwrap_or(&0) as f64;
    let max = *values.iter().max().unwrap_or(&0) as f64;
    max / (min + 1.0)
}

fn position_probabilities(base: char) -> [f64; 10] {
    match base {
        'A' => [0.51, 0.55, 0.57, 0.52, 0.48, 0.58, 0.57, 0.54, 0.50, 0.36],
        'C' => [0.29, 0.44, 0.55, 0.49, 0.52, 0.60, 0.60, 0.56, 0.51, 0.38],
        'G' => [0.62, 0.67, 0.74, 0.65, 0.61, 0.62, 0.52, 0.41, 0.31, 0.17],
        'T' => [0.51, 0.60, 0.69, 0.64, 0.62, 0.67, 0.58, 0.48, 0.39, 0.24],
        _ => [0.0; 10],
    }
}

fn content_probabilities(base: char) -> [f64; 10] {
    match base {
        'A' => [0.40, 0.55, 0.58, 0.58, 0.52, 0.48, 0.45, 0.45, 0.38, 0.19],
        'C' => [0.50, 0.63, 0.59, 0.50, 0.46, 0.45, 0.47, 0.56, 0.59, 0.33],
        'G' => [0.21, 0.40, 0.47, 0.50, 0.52, 0.56, 0.57, 0.52, 0.44, 0.23],
        'T' => [0.30, 0.49, 0.56, 0.53, 0.48, 0.48, 0.52, 0.57, 0.60, 0.51],
        _ => [0.0; 10],
    }
}

fn position_weight(base: char) -> f64 {
    match base {
        'A' => 0.062,
        'C' => 0.093,
        'G' => 0.205,
        'T' => 0.154,
        _ => 0.0,
    }
}

fn content_weight(base: char) -> f64 {
    match base {
        'A' => 0.084,
        'C' => 0.076,
        'G' => 0.081,
        'T' => 0.055,
        _ => 0.0,
    }
}

#[cfg(test)]
mod tests {
    use super::Fickett;

    #[test]
    fn gives_zero_for_tiny_sequences() {
        assert_eq!(Fickett::score("A"), 0.0);
    }
}
