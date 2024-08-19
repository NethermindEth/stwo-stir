//! This is the implementation of a circle modification of STIR algorithm from scratch,
//! moved from Python trying to keep it close to the original implementation.
use itertools::Itertools;
use num_bigint::{BigInt, Sign};
use num_traits::{One, Pow, ToPrimitive};
use crate::core::fields::cm31::CM31;
use super::*;
use super::fft::*;
use super::merkle_trees::*;
use super::poly_utils::*;

#[derive(Debug)]
struct Parameters<T> {
    root_of_unity: T,
    maxdeg_plus_1: u32,
    folding_params: Vec<usize>,
    eval_offsets: Vec<T>,
    eval_sizes: Vec<usize>,
    repetition_params: Vec<usize>,
    ood_rep: usize,
    modulus: u32,
    prim_root: T,
}

#[derive(Debug)]
struct Proof(Vec<u8>);

type GaussianF = CM31;

/// Generate an STIR proof that the polynomial that has the specified
/// values at successive powers of the specified root of unity has a
/// degree lower than maxdeg\_plus\_1
///
/// We use maxdeg\+1 instead of maxdeg because it's more mathematically
/// convenient in this case.
fn prove_low_degree(values: &[i128], params: &Parameters<GaussianF>, is_fake: bool) -> Proof {
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
    let mut xs = get_power_cycle(rt, params.eval_offsets[0]);
    {
        let xs_added =
            xs.iter().cloned().map(|x| x.conj(1)).collect_vec();
        xs.extend(&xs_added);
    }
    assert_eq!(values.len(), xs.len());
    assert_eq!(values.len(), params.eval_sizes[0] * 2);

    // Put the values into a Merkle tree. This is the root that the
    // proof will be checked against
    let mut m = merkelize(&vals);

    // Select a pseudo-random field element
    // r_fold = f.exp(prim_root, int.from_bytes(m[1], 'big') % (modulus + 1))
    let mut r_fold = {
        let pow = BigInt::from_bytes_be(Sign::Plus, &m[1]).rem_euclid(&BigInt::from(params.modulus + 1)).to_u32().unwrap();
        params.prim_root.pow(pow)
    };
    output.extend_from_slice(&m[1]);

    let mut last_g_hat: Option<Vec<i128>> = None;
    let mut last_folded_len: Option<usize> = None;

    for i in 1..=params.folding_params.len() {
        // fold using r-fold
        assert_eq!(params.eval_sizes[i - 1].rem_euclid(params.folding_params[i - 1]), 0);
        let folded_len = params.eval_sizes[i - 1] / params.folding_params[i - 1];
        last_folded_len = Some(folded_len);

        // should replace lagrange_interp with an fft?
        let rt2 = rt.pow(folded_len as u32);
        let xs2s = {
            let mut res = vec![get_power_cycle(rt2, GaussianF::one())];
            res.push(res[0].clone());
            (&mut res[1][1..]).reverse();
            res
        };

        let x_polys: Vec<Vec<Vec<i128>>> =
            (0..=1).flat_map(|l| {
                (0..folded_len).map(|k| {
                    let temp = f.circ_lagrange_interp(
                        &xs2s[l],
                        &(0..params.folding_params[i - 1])
                            .map(|j| vals[k + folded_len * j + params.eval_sizes[i - 1] * l])
                            .collect_vec(),
                        true,
                    );
                    temp
                }).collect_vec()
            })
                .collect_vec();

        let g_hat: Vec<i128> =
            (0..=1).flat_map(|l| {
                (0..folded_len).map(|k| {
                    let mul = r_fold * xs[k + params.eval_sizes[i - 1] * l].conj(1);
                    f.eval_circ_poly_at(&x_polys[k + folded_len * l], mul)
                }).collect_vec()
            })
                .collect_vec();
        last_g_hat = Some(g_hat.clone());

        if i == params.folding_params.len() {
            break;
        }

        assert_eq!(params.eval_sizes[i] % folded_len, 0);
        let expand_factor = params.eval_sizes[i] / folded_len;
        assert_eq!(params.eval_sizes[i - 1] % params.eval_sizes[i], 0);
        let eval_size_scale = params.eval_sizes[i - 1] / params.eval_sizes[i];

        rt = rt.pow(eval_size_scale as u32);
        let rt2 = rt.pow(expand_factor as u32);
        let p_offset = params.eval_offsets[i - 1].pow(params.folding_params[i - 1] as u32);

        let g_hat_shift = shift_domain(
            &g_hat,
            params.modulus,
            rt,
            p_offset,
            params.eval_offsets[i],
            expand_factor as u32,
        );
        let m2 = merkelize(&g_hat_shift);
        output.extend_from_slice(&m2[1]);

        xs = get_power_cycle(rt, params.eval_offsets[i]);
        xs.extend(xs.iter().map(|x| x.conj(1)).collect_vec());

        let r_outs = get_pseudorandom_element_outside_coset_circle(
            &output,
            params.modulus,
            params.prim_root,
            params.eval_sizes[i],
            params.eval_offsets[i],
            params.ood_rep,
        );

        let betas: Vec<i128> = r_outs
            .iter()
            .map(|&r| inv_fft_at_point(&g_hat, params.modulus, rt2, p_offset, r))
            .collect_vec();
        output.extend(betas.iter().flat_map(|&b| to_32_be_bytes(b)));

        r_fold = params.prim_root.pow(get_pseudorandom_indices(&output, params.modulus + 1, 1, 0, &mut vec![])[0] as u32);
        let r_comb = params.prim_root.pow(get_pseudorandom_indices(&output, params.modulus + 1, 1, 1, &mut vec![])[0] as u32);
        let t_vals = get_pseudorandom_indices(&output, 2 * folded_len as u32, params.repetition_params[i - 1], 2, &mut vec![]);
        let t_shifts = t_vals.iter().map(|&t| t / 2).collect_vec();
        let t_conj = t_vals.iter().map(|&t| t.rem_euclid(2)).collect_vec();
        assert_eq!((params.ood_rep + params.repetition_params[i - 1]).rem_euclid(2), 0);

        let rs: Vec<GaussianF> = r_outs.iter().cloned()
            .chain(t_shifts.iter().zip(t_conj.iter())
                .map(|(&t, &k)| (p_offset * rt2.pow(t as u32)).conj(k as u64)))
            .collect_vec();
        let g_rs: Vec<i128> = betas.iter().cloned()
            .chain(t_shifts.iter().zip(t_conj.iter())
                .map(|(&t, &k)| g_hat[t + k * folded_len]))
            .collect_vec();

        // maybe use compressed proof
        let oracle_branches =
            make_oracle_branches(params.folding_params[i - 1], params.eval_sizes[i - 1],
                                 &t_shifts, &t_conj, folded_len, &m);
        output.extend(oracle_branches);
        m = m2;

        let pol = f.circ_lagrange_interp(&rs, &g_rs, false);
        let pol_vals = xs.iter().map(|&x| f.eval_circ_poly_at(&pol, x)).collect_vec();
        let zpol = f.circ_zpoly(&rs, None);

        vals = (0..(2 * params.eval_sizes[i]))
            .map(|j| {
                f.div(
                    g_hat_shift[j] - pol_vals[j],
                    f.eval_circ_poly_at(&zpol, xs[j]),
                ) * f.geom_sum((xs[j] * r_comb).x(), rs.len() as u64)
            })
            .collect_vec();
    }

    let g_hat = last_g_hat.unwrap();
    let folded_len = last_folded_len.unwrap();
    let last_folding_param = params.folding_params.last().cloned().unwrap();
    let last_eval_offset = params.eval_offsets.last().cloned().unwrap();
    let g_pol = fft(&g_hat, params.modulus,
                    rt.pow(last_folding_param as u32),
                    last_eval_offset.pow(last_folding_param as u32));
    let final_deg = params.maxdeg_plus_1 as usize / params.folding_params.iter().product::<usize>();
    if !is_fake {
        assert!(g_pol[(2 * final_deg + 1)..].iter().all(|&x| x == 0));
    }
    output.extend(g_pol[..(2 * final_deg + 1)].iter().flat_map(|&c| to_32_be_bytes(c)));

    let t_vals = get_pseudorandom_indices(&output, 2 * folded_len as u32, params.repetition_params.last().cloned().unwrap(), 0, &mut vec![]);
    let t_shifts = t_vals.iter().map(|&t| t / 2).collect_vec();
    let t_conj = t_vals.iter().map(|&t| t.rem_euclid(2)).collect_vec();

    let oracle_branches =
        make_oracle_branches(last_folding_param,
                             params.eval_sizes.last().cloned().unwrap(),
                             &t_shifts, &t_conj, folded_len, &m);
    output.extend(oracle_branches);

    Proof(output)
}

fn make_oracle_branches(
    folding_param: usize,
    eval_size: usize,
    t_shifts: &[usize],
    t_conj: &[usize],
    folded_len: usize,
    m: &[Vec<u8>],
) -> Vec<u8> {
    let mut result = vec![];
    for (&t, &k) in t_shifts.iter().zip(t_conj.iter()) {
        for j in 0..folding_param {
            let branch = mk_branch(&m, t + j * folded_len + k * eval_size);
            result.extend(branch.into_iter().flatten())
        }
    }
    result
}

fn verify_low_degree_proof(proof: &Proof, params: &Parameters<GaussianF>) -> bool {
    macro_rules! reject_unless_eq {
        ($lhs:expr,$rhs:expr) => { if $lhs != $rhs { return false; } };
    }
    let f = PrimeField::new(params.modulus);

    reject_unless_eq!(params.folding_params.len(), params.eval_offsets.len());
    reject_unless_eq!(params.folding_params.len(), params.eval_sizes.len());
    reject_unless_eq!(params.folding_params.len(), params.repetition_params.len());

    let mut rt = params.root_of_unity;

    reject_unless_eq!(rt.pow(params.eval_sizes[0] as u32), GaussianF::one());
    reject_unless_eq!(rt.pow(params.eval_sizes[0] as u32 / 2), GaussianF::new(params.modulus as u64 - 1, 0));

    let mut proof_pos = 0;
    let m_root = &proof.0[proof_pos..(proof_pos + 32)];
    proof_pos += 32;
    let mut r_fold =
        params.prim_root.pow(
            BigInt::from_bytes_be(Sign::Plus, m_root)
                .rem_euclid(&BigInt::from(params.modulus + 1))
                .to_u32().unwrap(),
        );

    let mut pol: Option<Vec<Vec<i128>>> = None;
    let mut rs: Option<Vec<GaussianF>> = None;
    let mut zpol: Option<Vec<Vec<i128>>> = None;
    let mut r_comb: Option<GaussianF> = None;
    let mut m_root: Option<&[u8]> = None;

    for i in 1..params.folding_params.len() {
        reject_unless_eq!(params.eval_sizes[i - 1] % params.folding_params[i - 1], 0);
        let folded_len = params.eval_sizes[i - 1] / params.folding_params[i - 1];
        reject_unless_eq!(params.eval_sizes[i] % folded_len, 0);
        let expand_factor = params.eval_sizes[i] / folded_len;
        reject_unless_eq!(params.eval_sizes[i - 1] % params.eval_sizes[i], 0);
        let eval_size_scale = params.eval_sizes[i - 1] / params.eval_sizes[i];

        let rt_new = rt.pow(eval_size_scale as u32);
        let rt2 = rt_new.pow(expand_factor as u32);
        let p_offset = params.eval_offsets[i - 1].pow(params.folding_params[i - 1] as u32);

        let m2_root = &proof.0[proof_pos..proof_pos + 32];
        proof_pos += 32;

        let r_outs = get_pseudorandom_element_outside_coset_circle(
            &proof.0[..proof_pos], params.modulus, params.prim_root,
            params.eval_sizes[i], params.eval_offsets[i], params.ood_rep,
        );
        let betas: Vec<i128> = (0_usize..params.ood_rep).map(|j| {
            BigInt::from_bytes_be(Sign::Plus, &proof.0[(proof_pos + j * 32)..(proof_pos + (j + 1) * 32)])
                .rem_euclid(&BigInt::from(params.modulus))
                .to_i128().unwrap()
        }).collect();
        proof_pos += params.ood_rep * 32;

        let r_fold_new = params.prim_root.pow(get_pseudorandom_indices(&proof.0[..proof_pos], params.modulus + 1, 1, 0, &mut vec![])[0] as u32);
        let r_comb_new = params.prim_root.pow(get_pseudorandom_indices(&proof.0[..proof_pos], params.modulus + 1, 1, 1, &mut vec![])[0] as u32);

        let t_vals = get_pseudorandom_indices(&proof.0[..proof_pos], 2 * folded_len as u32, params.repetition_params[i - 1], 2, &mut vec![]);
        let t_shifts = t_vals.iter().map(|&t| t / 2).collect_vec();
        let t_conj = t_vals.iter().map(|&t| t.rem_euclid(2)).collect_vec();
        let rs_new: Vec<GaussianF> = r_outs.iter().cloned()
            .chain(t_shifts.iter().zip(t_conj.iter())
                .map(|(&t, &k)| (p_offset * rt2.pow(t as u32)).conj(k as u64)))
            .collect_vec();

        // should we change to storing the log since it should always be a power of 2?
        let branch_len = (params.eval_sizes[i - 1] as f64).log2().ceil() as usize + 2;
        let num_branches = params.folding_params[i - 1] * params.repetition_params[i - 1];
        let oracle_data: Vec<&[u8]> = (0_usize..(branch_len * num_branches)).map(|j| {
            &proof.0[(proof_pos + j * 32)..(proof_pos + (j + 1) * 32)]
        }).collect();
        proof_pos += branch_len * num_branches * 32;
        let oracle_branches: Vec<Vec<&[u8]>> = (0..num_branches).map(|j| {
            oracle_data[(j * branch_len)..((j + 1) * branch_len)].to_vec()
        }).collect();

        let rt2 = rt.pow(folded_len as u32);
        let xs2s = {
            let mut res = vec![get_power_cycle(rt2, GaussianF::new(1, 0))];
            res.push(res[0].clone());
            (&mut res[1][1..]).reverse();
            res
        };

        let mut g_hat = vec![];
        for k in 0_usize..params.repetition_params[i - 1] {
            let mut vals = vec![];
            let x0 = (rt.pow(t_shifts[k] as u32) * params.eval_offsets[i - 1])
                .conj(t_conj[k] as u64);
            for j in 0_usize..params.folding_params[i - 1] {
                let branch = &oracle_branches[k * params.folding_params[i - 1] + j];
                let ind = t_shifts[k] + j * folded_len;
                if let Some(m_root) = m_root {
                    verify_branch(m_root, ind + t_conj[k] * params.eval_sizes[i - 1], &branch);
                }
                let val = BigInt::from_bytes_be(Sign::Plus, branch[0])
                    .rem_euclid(&BigInt::from(params.modulus))
                    .to_i128().unwrap();

                let new_val = match pol {
                    Some(ref pol) => {
                        let rs = rs.as_ref().unwrap();
                        let r_comb = r_comb.unwrap();
                        let zpol = zpol.as_ref().unwrap();

                        let x = (rt.pow(ind as u32) * params.eval_offsets[i - 1])
                            .conj(t_conj[k] as u64);

                        f.div(val - f.eval_circ_poly_at(pol, x),
                              f.eval_circ_poly_at(zpol, x)) * f.geom_sum((x * r_comb).x(), rs.len() as u64)
                    }
                    None => val,
                };
                vals.push(new_val);
            }
            let new_g_hat = f.eval_circ_poly_at(&f.circ_lagrange_interp(&xs2s[t_conj[k]], &vals, true), r_fold * x0.conj(1));
            g_hat.push(new_g_hat);
        }

        rt = rt_new;
        m_root = Some(m2_root);
        r_fold = r_fold_new;
        r_comb = Some(r_comb_new);
        zpol = Some(f.circ_zpoly(&rs_new, None));
        let g_rs = betas.iter().cloned().chain(g_hat.iter().cloned()).collect_vec();
        pol = Some(f.circ_lagrange_interp(&rs_new, &g_rs, false));
        rs = Some(rs_new);
    }

    let final_deg = params.maxdeg_plus_1 as usize / params.folding_params.iter().product::<usize>();
    let g_pol: Vec<i128> = (0..(2 * final_deg + 1)).map(|j| {
        BigInt::from_bytes_be(Sign::Plus, &proof.0[(proof_pos + j * 32)..(proof_pos + (j + 1) * 32)])
            .rem_euclid(&BigInt::from(params.modulus))
            .to_i128().unwrap()
    }).collect();
    proof_pos += (2 * final_deg + 1) * 32;

    let last_eval_size = *params.eval_sizes.last().unwrap();
    let last_eval_offset = *params.eval_offsets.last().unwrap();
    let last_folding_param = *params.folding_params.last().unwrap();
    let last_repetition_param = *params.repetition_params.last().unwrap();

    reject_unless_eq!(last_eval_size % last_folding_param, 0);
    let folded_len = last_eval_size / last_folding_param;
    let rt2 = rt.pow(last_folding_param as u32);

    let t_vals = get_pseudorandom_indices(&proof.0[..proof_pos], 2 * folded_len as u32, params.repetition_params.last().cloned().unwrap(), 0, &mut vec![]);
    let t_shifts = t_vals.iter().map(|&t| t / 2).collect_vec();
    let t_conj = t_vals.iter().map(|&t| t.rem_euclid(2)).collect_vec();

    // Warning: code below copy-and-pasted from inside main loop
    let branch_len = (last_eval_size as f64).log2().ceil() as usize + 2;
    let num_branches = last_folding_param * last_repetition_param;
    let oracle_data: Vec<&[u8]> =
        (0..(branch_len * num_branches)).map(|j| &proof.0[(proof_pos + j * 32)..(proof_pos + (j + 1) * 32)]).collect();
    // proof_pos += branch_len * num_branches * 32;
    let oracle_branches: Vec<Vec<&[u8]>> = (0..num_branches).map(|j| {
        oracle_data[(j * branch_len)..((j + 1) * branch_len)].to_vec()
    }).collect();

    let rs = rs.as_ref().unwrap();

    let rt3 = rt.pow(folded_len as u32);
    let xs2s = {
        let mut res = vec![get_power_cycle(rt3, GaussianF::new(1, 0))];
        res.push(res[0].clone());
        (&mut res[1][1..]).reverse();
        res
    };

    for k in 0_usize..last_repetition_param {
        let mut vals = vec![];
        let x0 = (rt.pow(t_shifts[k] as u32) * last_eval_offset)
            .conj(t_conj[k] as u64);
        for j in 0..last_folding_param {
            let branch = &oracle_branches[k * last_folding_param + j];
            let ind = t_shifts[k] + j * folded_len;
            verify_branch(m_root.as_ref().unwrap(), ind + t_conj[k] * last_eval_size, branch);
            let val = BigInt::from_bytes_be(Sign::Plus, branch[0])
                .rem_euclid(&BigInt::from(params.modulus))
                .to_i128().unwrap();

            let x = (rt.pow(ind as u32) * last_eval_offset)
                .conj(t_conj[k] as u64);
            vals.push(f.div(
                val - f.eval_circ_poly_at(&pol.as_ref().unwrap(), x),
                f.eval_circ_poly_at(zpol.as_ref().unwrap(), x),
            ) * f.geom_sum((x * r_comb.unwrap()).x(), rs.len() as u64));
        }

        reject_unless_eq!(f.eval_circ_poly_at(&f.circ_lagrange_interp(&xs2s[t_conj[k]], &vals, true), r_fold * x0.conj(1)),
                          fft_inv(&g_pol, params.modulus, GaussianF::new(1,0),
                                  (rt2.pow(t_shifts[k] as u32) * last_eval_offset.pow(last_folding_param as u32))
                                      .conj(t_conj[k] as u64))[0]);
    }

    println!("STIR proof verified");
    true
}

fn get_pseudorandom_element_outside_coset_circle<F: KindaField>(
    seed: &[u8],
    modulus: u32,
    prim_root: F,
    coset_size: usize,
    shift: F,
    count: usize,
) -> Vec<F> {
    let adjust = 1;
    let m2 = modulus + 1;
    assert_eq!((modulus as i128 + adjust) % (coset_size as i128), 0);
    let cofactor = (modulus as i128 + adjust) / (coset_size as i128);
    let mut exclude = vec![];
    let mut ans = vec![];
    let mut start = 0;
    while ans.len() < count {
        let r = get_pseudorandom_indices(seed, m2 - coset_size as u32, count, start, &mut exclude)[0];
        start += 1;
        exclude.push(r);
        let t = r + 1 + r / (cofactor as usize - 1);

        let val = prim_root.pow(t as u32) * shift;
        if (val * shift).pow(coset_size as u32) != F::one() {
            ans.push(val);
        }
    }
    ans
}

#[cfg(test)]
mod tests {
    use crate::core::fields::m31;
    use super::super::fft::fft_inv;
    use super::*;

    #[test]
    fn test_circle_stir() {
        const MODULUS: u32 = m31::P;
        let prim_root = GaussianF::new(311014874, 1584694829);
        let root_of_unity = prim_root.pow((MODULUS + 1) / 2_u32.pow(10 + 2));
        let log_d = 10;
        let params = generate_parameters(MODULUS, prim_root, root_of_unity, prim_root, log_d, 128, 4, 1.0, 3);
        println!("{:?}", params);

        // Pure STIR tests
        let poly: Vec<i128> = (0..2_i128.pow(log_d + 1)).collect();
        let evaluations = fft_inv(&poly, MODULUS, root_of_unity, prim_root);
        let proof = prove_low_degree(&evaluations, &params, false);
        assert!(verify_low_degree_proof(&proof, &params));

        // let fakedata: Vec<i64> = evaluations.iter().enumerate().map(|(i, &x)| {
        //     if pow_mod(3, x as u32, 4096) > 400 { i as i64 } else { 39 }
        // }).collect();
        // let proof2 = prove_low_degree(&fakedata, &params, true);
        // println!("fake proof created");
        // match verify_low_degree_proof(&proof2, &params) {
        //     true => panic!("Fake data passed FRI"),
        //     false => println!("fake proof rejected"),
        // }
    }

    fn generate_parameters(
        modulus: u32,
        prim_root: GaussianF,
        root_of_unity: GaussianF,
        init_offset: GaussianF,
        log_d: u32,
        _security_param: u64,
        log_stopping_degree: u64,
        proximity_param: f64,
        log_folding_param: u32,
    ) -> Parameters<GaussianF> {
        let m = ((log_d as u64 - log_stopping_degree) / log_folding_param as u64) as usize;
        let size_l = get_power_cycle(root_of_unity, GaussianF::new(1, 0)).len();
        assert!(size_l.is_power_of_two());
        assert!(proximity_param > 0.0 && proximity_param <= 1.0);

        let mut eval_sizes = vec![size_l];
        let mut eval_offsets = vec![init_offset];
        let mut rt = root_of_unity;

        for _ in 1..=m {
            eval_sizes.push(eval_sizes.last().unwrap() / 2);
            eval_offsets.push(rt * init_offset);
            rt = rt.pow(2);
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
