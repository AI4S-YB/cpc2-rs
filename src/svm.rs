use std::sync::OnceLock;

use crate::error::Cpc2Error;

const RANGE_TEXT: &str = include_str!("../data/cpc2.range");
const MODEL_TEXT: &str = include_str!("../data/cpc2.model");

#[derive(Debug)]
pub struct Prediction {
    pub predicted_label: i32,
    pub coding_probability: f64,
}

#[derive(Debug)]
pub struct Predictor {
    scaler: FeatureScaler,
    model: BinarySvmModel,
}

#[derive(Debug)]
struct FeatureScaler {
    lower: f64,
    upper: f64,
    min: [f64; 4],
    max: [f64; 4],
}

#[derive(Debug)]
struct BinarySvmModel {
    gamma: f64,
    rho: f64,
    labels: [i32; 2],
    prob_a: f64,
    prob_b: f64,
    support_vectors: Vec<SupportVector>,
}

#[derive(Debug)]
struct SupportVector {
    coefficient: f64,
    features: [f64; 4],
}

pub fn embedded_predictor() -> &'static Predictor {
    static PREDICTOR: OnceLock<Predictor> = OnceLock::new();
    PREDICTOR.get_or_init(|| Predictor::from_embedded().expect("embedded CPC2 model must parse"))
}

impl Predictor {
    fn from_embedded() -> Result<Self, Cpc2Error> {
        Ok(Self {
            scaler: FeatureScaler::parse(RANGE_TEXT)?,
            model: BinarySvmModel::parse(MODEL_TEXT)?,
        })
    }

    pub fn predict(&self, raw_features: [f64; 4]) -> Prediction {
        let scaled = self.scaler.scale(raw_features);
        let decision = self.model.decision_value(scaled);
        let probability_for_label0 =
            sigmoid_predict(decision, self.model.prob_a, self.model.prob_b).clamp(1e-7, 1.0 - 1e-7);
        let probabilities = multiclass_probability(&[
            vec![0.0, probability_for_label0],
            vec![1.0 - probability_for_label0, 0.0],
        ]);

        let predicted_label = if probabilities[0] >= probabilities[1] {
            self.model.labels[0]
        } else {
            self.model.labels[1]
        };

        let coding_probability = if self.model.labels[0] == 1 {
            probabilities[0]
        } else {
            probabilities[1]
        };

        Prediction {
            predicted_label,
            coding_probability,
        }
    }
}

impl FeatureScaler {
    fn parse(input: &str) -> Result<Self, Cpc2Error> {
        let mut lines = input.lines().map(str::trim).filter(|line| !line.is_empty());
        let kind = lines
            .next()
            .ok_or_else(|| Cpc2Error::Parse("missing range header".to_string()))?;
        if kind != "x" {
            return Err(Cpc2Error::Parse(format!(
                "unsupported range header: {kind}"
            )));
        }

        let bounds = lines
            .next()
            .ok_or_else(|| Cpc2Error::Parse("missing scaler bounds".to_string()))?;
        let mut bound_parts = bounds.split_whitespace();
        let lower = parse_f64(bound_parts.next(), "range lower bound")?;
        let upper = parse_f64(bound_parts.next(), "range upper bound")?;

        let mut min = [0.0; 4];
        let mut max = [0.0; 4];
        let mut seen = [false; 4];

        for line in lines {
            let mut parts = line.split_whitespace();
            let index = parse_usize(parts.next(), "feature index")?;
            let minimum = parse_f64(parts.next(), "feature minimum")?;
            let maximum = parse_f64(parts.next(), "feature maximum")?;
            if !(1..=4).contains(&index) {
                continue;
            }
            min[index - 1] = minimum;
            max[index - 1] = maximum;
            seen[index - 1] = true;
        }

        if seen.iter().any(|value| !value) {
            return Err(Cpc2Error::Parse(
                "range file does not define all four CPC2 features".to_string(),
            ));
        }

        Ok(Self {
            lower,
            upper,
            min,
            max,
        })
    }

    fn scale(&self, features: [f64; 4]) -> [f64; 4] {
        let mut scaled = [0.0; 4];

        for idx in 0..4 {
            scaled[idx] = if self.max[idx] == self.min[idx] {
                0.0
            } else if features[idx] == self.min[idx] {
                self.lower
            } else if features[idx] == self.max[idx] {
                self.upper
            } else {
                self.lower
                    + (self.upper - self.lower) * (features[idx] - self.min[idx])
                        / (self.max[idx] - self.min[idx])
            };
        }

        scaled
    }
}

impl BinarySvmModel {
    fn parse(input: &str) -> Result<Self, Cpc2Error> {
        let mut gamma = None;
        let mut rho = None;
        let mut labels = None;
        let mut prob_a = None;
        let mut prob_b = None;
        let mut support_vectors = Vec::new();
        let mut in_support_vectors = false;

        for line in input.lines().map(str::trim).filter(|line| !line.is_empty()) {
            if in_support_vectors {
                support_vectors.push(parse_support_vector(line)?);
                continue;
            }

            if line == "SV" {
                in_support_vectors = true;
                continue;
            }

            let mut parts = line.split_whitespace();
            let key = parts
                .next()
                .ok_or_else(|| Cpc2Error::Parse("empty model line".to_string()))?;

            match key {
                "gamma" => gamma = Some(parse_f64(parts.next(), "gamma")?),
                "rho" => rho = Some(parse_f64(parts.next(), "rho")?),
                "label" => {
                    labels = Some([
                        parse_i32(parts.next(), "first model label")?,
                        parse_i32(parts.next(), "second model label")?,
                    ])
                }
                "probA" => prob_a = Some(parse_f64(parts.next(), "probA")?),
                "probB" => prob_b = Some(parse_f64(parts.next(), "probB")?),
                _ => {}
            }
        }

        Ok(Self {
            gamma: gamma.ok_or_else(|| Cpc2Error::Parse("model missing gamma".to_string()))?,
            rho: rho.ok_or_else(|| Cpc2Error::Parse("model missing rho".to_string()))?,
            labels: labels.ok_or_else(|| Cpc2Error::Parse("model missing labels".to_string()))?,
            prob_a: prob_a.ok_or_else(|| Cpc2Error::Parse("model missing probA".to_string()))?,
            prob_b: prob_b.ok_or_else(|| Cpc2Error::Parse("model missing probB".to_string()))?,
            support_vectors,
        })
    }

    fn decision_value(&self, scaled_features: [f64; 4]) -> f64 {
        self.support_vectors
            .iter()
            .map(|support_vector| {
                support_vector.coefficient
                    * rbf_kernel(self.gamma, support_vector.features, scaled_features)
            })
            .sum::<f64>()
            - self.rho
    }
}

fn parse_support_vector(line: &str) -> Result<SupportVector, Cpc2Error> {
    let mut parts = line.split_whitespace();
    let coefficient = parse_f64(parts.next(), "support vector coefficient")?;
    let mut features = [0.0; 4];

    for part in parts {
        let (raw_index, raw_value) = part
            .split_once(':')
            .ok_or_else(|| Cpc2Error::Parse(format!("invalid support vector feature: {part}")))?;
        let index = raw_index.parse::<usize>().map_err(|_| {
            Cpc2Error::Parse(format!("invalid support vector feature index: {raw_index}"))
        })?;
        let value = raw_value.parse::<f64>().map_err(|_| {
            Cpc2Error::Parse(format!("invalid support vector feature value: {raw_value}"))
        })?;
        if (1..=4).contains(&index) {
            features[index - 1] = value;
        }
    }

    Ok(SupportVector {
        coefficient,
        features,
    })
}

fn rbf_kernel(gamma: f64, left: [f64; 4], right: [f64; 4]) -> f64 {
    let mut distance = 0.0;
    for idx in 0..4 {
        let diff = left[idx] - right[idx];
        distance += diff * diff;
    }
    (-gamma * distance).exp()
}

fn sigmoid_predict(decision_value: f64, a: f64, b: f64) -> f64 {
    let f_apb = decision_value * a + b;
    if f_apb >= 0.0 {
        (-f_apb).exp() / (1.0 + (-f_apb).exp())
    } else {
        1.0 / (1.0 + f_apb.exp())
    }
}

fn multiclass_probability(pairwise_probabilities: &[Vec<f64>]) -> Vec<f64> {
    let class_count = pairwise_probabilities.len();
    let max_iter = usize::max(100, class_count);
    let eps = 0.005 / class_count as f64;

    let mut probabilities = vec![1.0 / class_count as f64; class_count];
    let mut q = vec![vec![0.0; class_count]; class_count];
    let mut qp = vec![0.0; class_count];

    for t in 0..class_count {
        for j in 0..t {
            q[t][t] += pairwise_probabilities[j][t] * pairwise_probabilities[j][t];
            q[t][j] = q[j][t];
        }
        for j in (t + 1)..class_count {
            q[t][t] += pairwise_probabilities[j][t] * pairwise_probabilities[j][t];
            q[t][j] = -pairwise_probabilities[j][t] * pairwise_probabilities[t][j];
        }
    }

    for _ in 0..max_iter {
        let mut p_q_p = 0.0;
        for t in 0..class_count {
            qp[t] = 0.0;
            for j in 0..class_count {
                qp[t] += q[t][j] * probabilities[j];
            }
            p_q_p += probabilities[t] * qp[t];
        }

        let max_error = qp
            .iter()
            .map(|value| (value - p_q_p).abs())
            .fold(0.0, f64::max);
        if max_error < eps {
            break;
        }

        for t in 0..class_count {
            let diff = (-qp[t] + p_q_p) / q[t][t];
            probabilities[t] += diff;
            p_q_p = (p_q_p + diff * (diff * q[t][t] + 2.0 * qp[t])) / (1.0 + diff) / (1.0 + diff);
            for j in 0..class_count {
                qp[j] = (qp[j] + diff * q[t][j]) / (1.0 + diff);
                probabilities[j] /= 1.0 + diff;
            }
        }
    }

    probabilities
}

fn parse_f64(value: Option<&str>, field: &str) -> Result<f64, Cpc2Error> {
    value
        .ok_or_else(|| Cpc2Error::Parse(format!("missing {field}")))?
        .parse::<f64>()
        .map_err(|_| Cpc2Error::Parse(format!("invalid {field}")))
}

fn parse_i32(value: Option<&str>, field: &str) -> Result<i32, Cpc2Error> {
    value
        .ok_or_else(|| Cpc2Error::Parse(format!("missing {field}")))?
        .parse::<i32>()
        .map_err(|_| Cpc2Error::Parse(format!("invalid {field}")))
}

fn parse_usize(value: Option<&str>, field: &str) -> Result<usize, Cpc2Error> {
    value
        .ok_or_else(|| Cpc2Error::Parse(format!("missing {field}")))?
        .parse::<usize>()
        .map_err(|_| Cpc2Error::Parse(format!("invalid {field}")))
}

#[cfg(test)]
mod tests {
    use super::embedded_predictor;

    #[test]
    fn parses_embedded_model() {
        let predictor = embedded_predictor();
        let prediction = predictor.predict([100.0, 0.4, 8.0, 1.0]);
        assert!((0.0..=1.0).contains(&prediction.coding_probability));
    }
}
