use num_traits::Zero;

use crate::{identity, Matrix};

pub fn supplement_basis(m: Matrix) -> Option<Matrix> {
    let n = m.len();
    let k = m.first().map(|row| row.len()).unwrap_or(0);
    if k > n || m.iter().any(|row| row.len() != k) {
        return None;
    }

    let mut m_work = m.clone();
    let mut b = identity(n);

    for s in 0..k {
        let mut t_opt = None;
        for t in s..n {
            if !m_work[t][s].is_zero() {
                t_opt = Some(t);
                break;
            }
        }
        let t = match t_opt {
            Some(idx) => idx,
            None => return None,
        };
        let d = m_work[t][s].clone().recip();

        if t != s {
            for row in 0..n {
                b[row][t] = b[row][s].clone();
            }
        }
        for row in 0..n {
            b[row][s] = m[row][s].clone();
        }

        for j in (s + 1)..k {
            let msj = &d * &m_work[t][j];
            m_work[s][j] = msj.clone();
            if t != s {
                m_work[t][j] = msj.clone();
            }
            for i in 0..n {
                if i == s || i == t {
                    continue;
                }
                let adjustment = m_work[i][s].clone() * m_work[t][j].clone();
                m_work[i][j] = m_work[i][j].clone() - adjustment;
            }
        }
    }

    Some(b)
}
