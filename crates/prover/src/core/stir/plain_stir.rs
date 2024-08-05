//! This is the implementation of an original STIR (non-circle) algorithm from scratch,
//! moved from Pyhton trying to keep it close to the original implementation.

use num_bigint::{BigInt, Sign};
use num_traits::ToPrimitive;
use super::*;
use super::fft::*;
use super::merkle_trees::*;
use super::poly_utils::*;

#[allow(dead_code)]
#[derive(Debug)]
struct Parameters {
    root_of_unity: i64,
    maxdeg_plus_1: u32,
    folding_params: Vec<usize>,
    eval_offsets: Vec<i64>,
    eval_sizes: Vec<usize>,
    repetition_params: Vec<usize>,
    ood_rep: usize,
    modulus: u32,
    prim_root: u32,
}

#[derive(Debug)]
struct Proof(Vec<u8>);

/// Generate an STIR proof that the polynomial that has the specified
/// values at successive powers of the specified root of unity has a
/// degree lower than maxdeg\_plus\_1
///
/// We use maxdeg\+1 instead of maxdeg because it's more mathematically
/// convenient in this case.
fn prove_low_degree(values: &[i64], params: &Parameters, is_fake: bool) -> Proof {
    let f = PrimeField::new(params.modulus);
    if is_fake {
        println!("Faking proof {} values are degree <= {}", values.len(), params.maxdeg_plus_1);
    } else {
        println!("Proving {} values are degree <= {}", values.len(), params.maxdeg_plus_1);
    }

    assert_eq!(params.folding_params.len(), params.eval_offsets.len());
    assert_eq!(params.folding_params.len(), params.eval_sizes.len());
    assert_eq!(params.folding_params.len(), params.repetition_params.len());

    let mut output = vec![];

    let mut rt = params.root_of_unity;
    let mut vals = values.to_vec();
    // Calculate the set of x coordinates
    let mut xs = get_power_cycle(rt, params.modulus, params.eval_offsets[0]);
    assert_eq!(values.len(), xs.len());
    assert_eq!(values.len(), params.eval_sizes[0] as usize);

    // Put the values into a Merkle tree. This is the root that the
    // proof will be checked against
    let mut m = merkelize(&vals);

    // Select a pseudo-random field element
    let mut r_fold = BigInt::from_bytes_be(Sign::Plus, &m[1]).rem_euclid(&BigInt::from(params.modulus)).to_i64().unwrap();

    output.extend_from_slice(&m[1]);

    let mut last_g_hat: Option<Vec<i64>> = None;
    let mut last_folded_len: Option<usize> = None;

    for i in 1..=params.folding_params.len() {
        // fold using r-fold
        assert_eq!(params.eval_sizes[i - 1].rem_euclid(params.folding_params[i - 1]), 0);
        let folded_len = params.eval_sizes[i - 1] / params.folding_params[i - 1];
        last_folded_len = Some(folded_len);
        // should replace lagrange_interp with an fft?
        let x_polys: Vec<Vec<i64>> = (0..folded_len)
            .map(|k| {
                f.lagrange_interp(
                    &(0..params.folding_params[i - 1])
                        .map(|j| xs[k + folded_len * j])
                        .collect::<Vec<_>>(),
                    &(0..params.folding_params[i - 1])
                        .map(|j| vals[k + folded_len * j])
                        .collect::<Vec<_>>(),
                )
            })
            .collect();
        let g_hat: Vec<i64> = x_polys.iter().map(|p| f.eval_poly_at(p, r_fold)).collect();
        last_g_hat = Some(g_hat.clone());

        if i == params.folding_params.len() {
            break;
        }

        assert_eq!(params.eval_sizes[i] % folded_len, 0);
        let expand_factor = params.eval_sizes[i] / folded_len;
        assert_eq!(params.eval_sizes[i - 1] % params.eval_sizes[i], 0);
        let eval_size_scale = params.eval_sizes[i - 1] / params.eval_sizes[i];

        rt = f.exp(rt, eval_size_scale as u32);
        let rt2 = f.exp(rt, expand_factor as u32);
        let p_offset = f.exp(params.eval_offsets[i - 1], params.folding_params[i - 1] as u32);
        let inv_offset = f.inv(p_offset as i64);
        let g_hat_shift = shift_domain(
            &g_hat,
            params.modulus,
            rt,
            f.mul(params.eval_offsets[i], inv_offset),
            expand_factor as u32,
        );
        let m2 = merkelize(&g_hat_shift);
        output.extend_from_slice(&m2[1]);

        xs = get_power_cycle(rt, params.modulus, params.eval_offsets[i]);

        let r_outs = get_pseudorandom_element_outside_coset(
            &output,
            params.modulus,
            params.prim_root as i64,
            params.eval_sizes[i],
            params.eval_offsets[i],
            params.ood_rep,
            false,
        );

        let betas: Vec<i64> = r_outs
            .iter()
            .map(|&r| inv_fft_at_point(&g_hat, params.modulus, rt2, f.mul(r, inv_offset)))
            .collect();
        output.extend(betas.iter().flat_map(|&b| to_32_be_bytes(b)));

        r_fold = get_pseudorandom_indices(&output, params.modulus, 1, 0)[0] as i64;
        // DegCor doesn't work well if the parameter is 0
        let r_comb = get_pseudorandom_indices(&output, params.modulus - 1, 1, 1)[0] + 1;
        let t_shifts = get_pseudorandom_indices(&output, folded_len as u32, params.repetition_params[i - 1], 2);

        let rs: Vec<i64> = r_outs
            .iter()
            .cloned()
            .chain(t_shifts.iter().map(|&t| f.mul(p_offset, f.exp(rt2, t as u32))))
            .collect();
        let g_rs: Vec<i64> = betas
            .iter()
            .cloned()
            .chain(t_shifts.iter().map(|&t| g_hat[t]))
            .collect();

        // maybe use compressed proof
        let oracle_branches =
            make_oracle_branches(&t_shifts, params.folding_params[i - 1], folded_len, &m);
        output.extend(oracle_branches);
        m = m2;

        let pol = f.lagrange_interp(&rs, &g_rs);
        let pol_vals = fft(&inv_shift_poly(&pol, params.modulus, params.eval_offsets[i]), params.modulus, rt, false);
        assert_eq!(
            pol_vals,
            xs.iter().map(|&x| f.eval_poly_at(&pol, x)).collect::<Vec<_>>()
        );

        vals = (0..params.eval_sizes[i])
            .map(|j| {
                f.mul(
                    f.div(
                        g_hat_shift[j] - pol_vals[j],
                        f.prod(&rs.iter().map(|&r| xs[j] - r).collect::<Vec<_>>()),
                    ),
                    f.geom_sum(xs[j] * r_comb as i64, rs.len() as u64),
                )
            })
            .collect();
    }

    let g_hat = last_g_hat.unwrap();
    let folded_len = last_folded_len.unwrap();
    let g_pol = fft(&g_hat, params.modulus, f.exp(rt, params.folding_params.last().cloned().unwrap() as u32), true);
    let final_deg = params.maxdeg_plus_1 as usize / params.folding_params.iter().product::<usize>();
    if !is_fake {
        assert!(g_pol[final_deg..].iter().all(|&x| x == 0));
    }
    output.extend(g_pol[..final_deg].iter().flat_map(|&c| to_32_be_bytes(c)));

    let t_shifts =
        get_pseudorandom_indices(&output, folded_len as u32, params.repetition_params.last().cloned().unwrap(), 0);
    let oracle_branches =
        make_oracle_branches(&t_shifts, params.folding_params.last().cloned().unwrap(), folded_len, &m);
    output.extend(oracle_branches);

    Proof(output)
}

fn make_oracle_branches(t_shifts: &[usize], last_folding_param: usize, folded_len: usize, m: &[Vec<u8>]) -> Vec<u8> {
    let mut result = vec![];
    for t in t_shifts {
        for j in 0..last_folding_param {
            let branch = mk_branch(&m, t + j * folded_len);
            result.extend(branch.into_iter().flatten())
        }
    }
    result
}

fn verify_low_degree_proof(proof: &Proof, params: &Parameters) -> bool {
    macro_rules! reject_unless_eq {
        ($lhs:expr,$rhs:expr) => { if $lhs != $rhs { return false; } };
    }
    let f = PrimeField::new(params.modulus);

    reject_unless_eq!(params.folding_params.len(), params.eval_offsets.len());
    reject_unless_eq!(params.folding_params.len(), params.eval_sizes.len());
    reject_unless_eq!(params.folding_params.len(), params.repetition_params.len());

    let mut rt = params.root_of_unity;

    reject_unless_eq!(f.exp(rt, params.eval_sizes[0] as u32), 1);
    reject_unless_eq!(f.exp(rt, params.eval_sizes[0] as u32 / 2), params.modulus as i64 - 1);

    let mut proof_pos = 0;
    let m_root = &proof.0[proof_pos..(proof_pos + 32)];
    proof_pos += 32;
    let mut r_fold = BigInt::from_bytes_be(Sign::Plus, m_root)
        .rem_euclid(&BigInt::from(params.modulus))
        .to_i64().unwrap();

    let mut pol: Option<Vec<i64>> = None;
    let mut rs: Option<Vec<i64>> = None;
    let mut r_comb: Option<i64> = None;
    let mut m_root: Option<&[u8]> = None;

    for i in 1..params.folding_params.len() {
        reject_unless_eq!(params.eval_sizes[i - 1] % params.folding_params[i - 1], 0);
        let folded_len = params.eval_sizes[i - 1] / params.folding_params[i - 1];
        reject_unless_eq!(params.eval_sizes[i] % folded_len, 0);
        let expand_factor = params.eval_sizes[i] / folded_len;
        reject_unless_eq!(params.eval_sizes[i - 1] % params.eval_sizes[i], 0);
        let eval_size_scale = params.eval_sizes[i - 1] / params.eval_sizes[i];

        let rt_new = f.exp(rt, eval_size_scale as u32);
        let rt2 = f.exp(rt_new, expand_factor as u32);
        let p_offset = f.exp(params.eval_offsets[i - 1], params.folding_params[i - 1] as u32);

        let m2_root = &proof.0[proof_pos..proof_pos + 32];
        proof_pos += 32;

        let r_outs = get_pseudorandom_element_outside_coset(
            &proof.0[..proof_pos], params.modulus, params.prim_root as i64,
            params.eval_sizes[i], params.eval_offsets[i], params.ood_rep,
            false,
        );
        let betas: Vec<i64> = (0_usize..params.ood_rep).map(|j| {
            BigInt::from_bytes_be(Sign::Plus, &proof.0[(proof_pos + j * 32)..(proof_pos + (j + 1) * 32)])
                .rem_euclid(&BigInt::from(params.modulus))
                .to_i64().unwrap()
        }).collect();
        proof_pos += params.ood_rep * 32;

        let r_fold_new = get_pseudorandom_indices(&proof.0[..proof_pos], params.modulus, 1, 0)[0] as i64;
        let r_comb_new = get_pseudorandom_indices(&proof.0[..proof_pos], params.modulus - 1, 1, 1)[0] as i64 + 1;
        let t_shifts = get_pseudorandom_indices(&proof.0[..proof_pos], folded_len as u32, params.repetition_params[i - 1], 2);
        let rs_new: Vec<i64> = r_outs.iter().cloned()
            .chain(t_shifts.iter().map(|&t| f.mul(p_offset, f.exp(rt2, t as u32))))
            .collect();

        // should we change to storing the log since it should always be a power of 2?
        let branch_len = (params.eval_sizes[i - 1] as f64).log2().ceil() as usize + 1;
        let num_branches = params.folding_params[i - 1] * params.repetition_params[i - 1];
        let oracle_data: Vec<&[u8]> = (0_usize..(branch_len * num_branches)).map(|j| {
            &proof.0[(proof_pos + j * 32)..(proof_pos + (j + 1) * 32)]
        }).collect();
        proof_pos += branch_len * num_branches * 32;
        let oracle_branches: Vec<Vec<&[u8]>> = (0..num_branches).map(|j| {
            oracle_data[(j * branch_len)..((j + 1) * branch_len)].to_vec()
        }).collect();

        let mut g_hat = vec![];
        for k in 0_usize..params.repetition_params[i - 1] {
            let mut vals = vec![];
            let mut xs = vec![];
            for j in 0_usize..params.folding_params[i - 1] {
                let branch = &oracle_branches[k * params.folding_params[i - 1] + j];
                let ind = t_shifts[k] + j * folded_len;
                if let Some(m_root) = m_root {
                    verify_branch(m_root, ind, &branch);
                }
                let val = BigInt::from_bytes_be(Sign::Plus, branch[0])
                    .rem_euclid(&BigInt::from(params.modulus))
                    .to_i64().unwrap();

                let x = f.mul(f.exp(rt, ind as u32), params.eval_offsets[i - 1]);
                xs.push(x);
                let new_val = match pol {
                    Some(ref pol) => {
                        let rs = rs.as_ref().unwrap();
                        let r_comb = r_comb.unwrap();
                        f.mul(f.div(val - f.eval_poly_at(&pol, x),
                                    f.prod(&rs.iter().map(|&r| x - r).collect::<Vec<_>>())),
                              f.geom_sum(x * r_comb, rs.len() as u64))
                    }
                    None => val,
                };
                vals.push(new_val);
            }
            g_hat.push(f.eval_poly_at(&f.lagrange_interp(&xs, &vals), r_fold));
        }

        rt = rt_new;
        m_root = Some(m2_root);
        r_fold = r_fold_new;
        r_comb = Some(r_comb_new);
        rs = Some(rs_new);
        let g_rs = betas.iter().cloned().chain(g_hat.iter().cloned()).collect::<Vec<_>>();
        pol = Some(f.lagrange_interp(rs.as_ref().unwrap(), &g_rs));
    }

    let final_deg = params.maxdeg_plus_1 as usize / params.folding_params.iter().product::<usize>();
    let g_pol: Vec<i64> = (0..final_deg).map(|j| {
        BigInt::from_bytes_be(Sign::Plus, &proof.0[(proof_pos + j * 32)..(proof_pos + (j + 1) * 32)])
            .rem_euclid(&BigInt::from(params.modulus))
            .to_i64().unwrap()
    }).collect();
    proof_pos += final_deg * 32;

    let last_eval_size = *params.eval_sizes.last().unwrap();
    let last_eval_offset = *params.eval_offsets.last().unwrap();
    let last_folding_param = *params.folding_params.last().unwrap();
    let last_repetition_param = *params.repetition_params.last().unwrap();

    reject_unless_eq!(last_eval_size % last_folding_param, 0);
    let folded_len = last_eval_size / last_folding_param;
    let rt2 = f.exp(rt, last_folding_param as u32);
    let t_shifts = get_pseudorandom_indices(&proof.0[..proof_pos], folded_len as u32, last_repetition_param, 0);

    // Warning: code below copy-and-pasted from inside main loop
    let branch_len = (last_eval_size as f64).log2().ceil() as usize + 1;
    let num_branches = last_folding_param * last_repetition_param;
    let oracle_data: Vec<&[u8]> =
        (0..(branch_len * num_branches)).map(|j| &proof.0[(proof_pos + j * 32)..(proof_pos + (j + 1) * 32)]).collect();
    // proof_pos += branch_len * num_branches * 32;
    let oracle_branches: Vec<Vec<&[u8]>> = (0..num_branches).map(|j| {
        oracle_data[(j * branch_len)..((j + 1) * branch_len)].to_vec()
    }).collect();

    let rs = rs.as_ref().unwrap();
    for k in 0_usize..last_repetition_param {
        let mut vals = vec![];
        let mut xs = vec![];
        for j in 0..last_folding_param {
            let branch = &oracle_branches[k * last_folding_param + j];
            let ind = t_shifts[k] + j * folded_len;
            verify_branch(m_root.as_ref().unwrap(), ind, branch);
            let val = BigInt::from_bytes_be(Sign::Plus, branch[0])
                .rem_euclid(&BigInt::from(params.modulus))
                .to_i64().unwrap();

            let x = f.mul(f.exp(rt, ind as u32), last_eval_offset);
            xs.push(x);
            vals.push(f.mul(f.div(val - f.eval_poly_at(&pol.as_ref().unwrap(), x),
                                  f.prod(&rs.iter().map(|&r| x - r).collect::<Vec<_>>())),
                            f.geom_sum(x * r_comb.unwrap(), rs.len() as u64)));
        }
        reject_unless_eq!(f.eval_poly_at(&f.lagrange_interp(&xs, &vals), r_fold),
                          f.eval_poly_at(&g_pol, f.exp(rt2, t_shifts[k] as u32)));
    }

    println!("STIR proof verified");
    true
}

fn get_pseudorandom_element_outside_coset(
    seed: &[u8],
    modulus: u32,
    prim_root: i64,
    coset_size: usize,
    shift: i64,
    count: usize,
    is_circle: bool /* default = false*/,
) -> Vec<i64> {
    let adjust = if is_circle { 1 } else { -1 };
    let m2 = if is_circle { modulus + 1 } else { modulus };
    assert_eq!((modulus as i64 + adjust) % (coset_size as i64), 0);
    let cofactor = (modulus as i64 + adjust) / (coset_size as i64);
    let rands = get_pseudorandom_indices(seed, m2 - coset_size as u32, count, 0);
    let mut ans = Vec::new();
    for r in rands {
        let t = r + 1 + r / (cofactor as usize - 1);
        if !is_circle && t == modulus as usize {
            ans.push(0);
        } else {
            ans.push((pow_mod(prim_root, t as u32, modulus) * shift).rem_euclid(modulus as i64));
        }
    }
    ans
}

fn get_power_cycle(r: i64, modulus: u32, offset: i64) -> Vec<i64> {
    let mut o = vec![offset];
    loop {
        let next = (o[o.len() - 1] * r) % (modulus as i64);
        if next == offset {
            break;
        }
        o.push(next);
    }
    o
}

#[cfg(test)]
mod tests {
    use super::super::fft::fft;
    use super::super::pow_mod;
    use super::*;

    #[test]
    fn test_stir_basics() {
        let modulus = 2_u32.pow(20) * 7 + 1;
        let prim_root = 3_u32;
        let root_of_unity = pow_mod(prim_root as i64, (modulus - 1) / 2_u32.pow(10 + 2), modulus);
        let log_d = 10;
        let params = generate_parameters(modulus, prim_root, root_of_unity, 1, log_d, 128, 4, 1.0, 3);
        println!("{:?}", params);

        // Pure STIR tests
        let poly: Vec<i64> = (0..2_i64.pow(log_d)).collect();
        let evaluations = fft(&poly, modulus, root_of_unity, false);
        let proof = prove_low_degree(&evaluations, &params, false);
        assert!(verify_low_degree_proof(&proof, &params));

        let fakedata: Vec<i64> = evaluations.iter().enumerate().map(|(i, &x)| {
            if pow_mod(3, x as u32, 4096) > 400 { i as i64 } else { 39 }
        }).collect();
        let proof2 = prove_low_degree(&fakedata, &params, true);
        println!("fake proof created");
        match verify_low_degree_proof(&proof2, &params) {
            true => panic!("Fake data passed FRI"),
            false => println!("fake proof rejected"),
        }
    }

    fn generate_parameters(
        modulus: u32,
        prim_root: u32,
        root_of_unity: i64,
        init_offset: u64,
        log_d: u32,
        _security_param: u64,
        log_stopping_degree: u64,
        proximity_param: f64,
        log_folding_param: u32,
    ) -> Parameters {
        let m = ((log_d as u64 - log_stopping_degree) / log_folding_param as u64) as usize;
        let size_l = get_power_cycle(root_of_unity, modulus, 1).len();
        assert!(size_l.is_power_of_two());
        assert!(proximity_param > 0.0 && proximity_param <= 1.0);

        let mut eval_sizes = vec![size_l as usize];
        let mut eval_offsets = vec![init_offset as i64];
        let mut rt = root_of_unity;

        for _ in 1..=m {
            eval_sizes.push(eval_sizes.last().unwrap() / 2);
            eval_offsets.push(rt);
            rt = pow_mod(rt, 2, modulus);
        }

        Parameters {
            root_of_unity,
            maxdeg_plus_1: 2_u32.pow(log_d),
            folding_params: vec![2_usize.pow(log_folding_param); m + 1],
            eval_offsets,
            eval_sizes,
            repetition_params: vec![3; m + 1],
            ood_rep: 1,
            modulus,
            prim_root,
        }
    }
}
