//! This is the implementation of a circle modification of STIR algorithm from scratch,
//! moved from Python trying to keep it close to the original implementation.

use itertools::Itertools;
use num_bigint::BigUint;
use num_traits::{ToPrimitive, Zero};
use crate::core::fields::{ComplexConjugate};
use super::*;
use super::fft;
use super::merkle_trees::*;
use super::poly_utils;

#[derive(Debug)]
struct Parameters<T> {
    root_of_unity: T,
    maxdeg_plus_1: u32,
    folding_params: Vec<usize>,
    eval_offsets: Vec<T>,
    eval_sizes: Vec<usize>,
    repetition_params: Vec<usize>,
    ood_rep: usize,
    prim_root: T,
}

#[derive(Debug)]
struct StirProofLayer {
    merkle_root: Hash,
    betas: Vec<M31>,
    oracle_branches: Vec<Vec<Hash>>,
}

#[derive(Debug)]
struct StirProof {
    merkle_root: Hash,
    layers: Vec<StirProofLayer>,
    last_g_pol: Vec<M31>,
    last_oracle_branches: Vec<Vec<Hash>>,

    /// Contains all proof bytes linearized, used as a Fiat-Shamir seed
    // TODO(alex): Remove once we switch from a makeshift hasher to a Blake2sChannel
    legacy_proof: Vec<u8>,
}

/// Generate an STIR proof that the polynomial that has the specified
/// values at successive powers of the specified root of unity has a
/// degree lower than maxdeg\_plus\_1
///
/// We use maxdeg\+1 instead of maxdeg because it's more mathematically
/// convenient in this case.
fn prove_low_degree<F: StirField<M31>>(values: &[M31], params: &Parameters<F>, is_fake: bool) -> StirProof {
    if is_fake {
        println!("Faking proof {} values are degree <= {}", values.len(), params.maxdeg_plus_1);
    } else {
        println!("Proving {} values are degree <= {}", values.len(), params.maxdeg_plus_1);
    }

    assert_eq!(params.folding_params.len(), params.eval_offsets.len());
    assert_eq!(params.folding_params.len(), params.eval_sizes.len());
    assert_eq!(params.folding_params.len(), params.repetition_params.len());

    let mut proof = StirProof {
        merkle_root: [0; 32],
        layers: vec![],
        last_g_pol: vec![],
        last_oracle_branches: vec![],
        legacy_proof: vec![],
    };

    let mut rt = params.root_of_unity;
    let mut vals = values.to_vec();
    // Calculate the set of x coordinates
    let mut xs = get_power_cycle(rt, params.eval_offsets[0]);
    {
        let xs_added =
            xs.iter().cloned().map(|x| x.complex_conjugate()).collect_vec();
        xs.extend(&xs_added);
    }
    assert_eq!(values.len(), xs.len());
    assert_eq!(values.len(), params.eval_sizes[0] * 2);

    // Put the values into a Merkle tree. This is the root that the
    // proof will be checked against
    let mut m = merkelize(&vals);

    // Select a pseudo-random field element
    let mut r_fold = {
        let pow = (BigUint::from_bytes_be(&m[1]) % BigUint::from(MODULUS + 1)).to_u128().unwrap();
        params.prim_root.pow(pow)
    };
    proof.legacy_proof.extend_from_slice(&m[1]);
    proof.merkle_root = m[1].clone().try_into().expect("Incorrect merkle root length");

    let mut last_g_hat: Option<Vec<M31>> = None;
    let mut last_folded_len: Option<usize> = None;

    for i in 1..=params.folding_params.len() {
        // fold using r-fold
        assert_eq!(params.eval_sizes[i - 1] % params.folding_params[i - 1], 0);
        let folded_len = params.eval_sizes[i - 1] / params.folding_params[i - 1];
        last_folded_len = Some(folded_len);

        // should replace lagrange_interp with an fft?
        let rt2 = rt.pow(folded_len as u128);
        let xs2s = {
            let mut res = vec![get_power_cycle(rt2, F::one())];
            res.push(res[0].clone());
            (&mut res[1][1..]).reverse();
            res
        };

        let x_polys: Vec<(Vec<M31>, Vec<M31>)> =
            (0..=1).flat_map(|l| {
                (0..folded_len).map(|k| {
                    let temp = poly_utils::circ_lagrange_interp(
                        &xs2s[l],
                        &(0..params.folding_params[i - 1])
                            .map(|j| vals[k + folded_len * j + params.eval_sizes[i - 1] * l])
                            .collect_vec(),
                        true,
                    );
                    temp
                }).collect_vec()
            }).collect_vec();

        let g_hat: Vec<M31> =
            (0..=1).flat_map(|l| {
                (0..folded_len).map(|k| {
                    let mul = r_fold * xs[k + params.eval_sizes[i - 1] * l].complex_conjugate();
                    poly_utils::eval_circ_poly_at(&x_polys[k + folded_len * l], &mul)
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

        rt = rt.pow(eval_size_scale as u128);
        let rt2 = rt.pow(expand_factor as u128);
        let p_offset = params.eval_offsets[i - 1].pow(params.folding_params[i - 1] as u128);

        let g_hat_shift = fft::shift_domain(
            &g_hat,
            rt,
            p_offset,
            params.eval_offsets[i],
            expand_factor as u32,
        );
        let m2 = merkelize(&g_hat_shift);
        proof.legacy_proof.extend_from_slice(&m2[1]);

        xs = get_power_cycle(rt, params.eval_offsets[i]);
        xs.extend(xs.iter().map(|x| x.complex_conjugate()).collect_vec());

        let r_outs = get_pseudorandom_element_outside_coset_circle(
            &proof.legacy_proof,
            params.prim_root,
            params.eval_sizes[i],
            params.eval_offsets[i],
            params.ood_rep,
        );

        let betas: Vec<M31> = r_outs
            .iter()
            .map(|&r| fft::inv_fft_at_point(&g_hat, rt2, p_offset, r))
            .collect_vec();
        proof.legacy_proof.extend(betas.iter().flat_map(|b| to_32_be_bytes(b.0 as u128)));

        r_fold = params.prim_root.pow(get_pseudorandom_indices(&proof.legacy_proof, MODULUS + 1, 1, 0, &mut vec![])[0] as u128);
        let r_comb = params.prim_root.pow(get_pseudorandom_indices(&proof.legacy_proof, MODULUS + 1, 1, 1, &mut vec![])[0] as u128);
        let t_vals = get_pseudorandom_indices(&proof.legacy_proof, 2 * folded_len as u32, params.repetition_params[i - 1], 2, &mut vec![]);
        let t_shifts = t_vals.iter().map(|&t| t / 2).collect_vec();
        let t_conj = t_vals.iter().map(|&t| t % 2).collect_vec();
        assert_eq!((params.ood_rep + params.repetition_params[i - 1]) % 2, 0);

        let rs: Vec<F> = r_outs.iter().cloned()
            .chain(t_shifts.iter().zip(t_conj.iter())
                .map(|(&t, &k)| {
                    conjugate_with_parity(p_offset * rt2.pow(t as u128), k)
                }))
            .collect_vec();

        // maybe use compressed proof
        let oracle_branches =
            make_oracle_branches(params.folding_params[i - 1], params.eval_sizes[i - 1],
                                 &t_shifts, &t_conj, folded_len, &m);
        proof.legacy_proof.extend(oracle_branches.iter().cloned().flatten().flatten().collect_vec());

        let proof_layer = StirProofLayer {
            merkle_root: m2[1].clone().try_into().expect("Incorrect merkle root length"),
            betas: betas.clone(),
            oracle_branches,
        };
        proof.layers.push(proof_layer);

        let g_rs: Vec<M31> = betas.into_iter()
            .chain(t_shifts.iter().zip(t_conj.iter())
                .map(|(&t, &k)| g_hat[t + k * folded_len]))
            .collect_vec();
        let pol = poly_utils::circ_lagrange_interp(&rs, &g_rs, false);
        let pol_vals = xs.iter().map(|&x| poly_utils::eval_circ_poly_at(&pol, &x)).collect_vec();
        let zpol = poly_utils::circ_zpoly(&rs, None);

        vals = (0..(2 * params.eval_sizes[i]))
            .map(|j| {
                (g_hat_shift[j] - pol_vals[j]) / poly_utils::eval_circ_poly_at(&zpol, &xs[j]) * poly_utils::geom_sum((xs[j] * r_comb).x(), rs.len() as u64)
            })
            .collect_vec();

        m = m2;
    }

    let g_hat = last_g_hat.unwrap();
    let folded_len = last_folded_len.unwrap();
    let last_folding_param = params.folding_params.last().cloned().unwrap();
    let last_eval_offset = params.eval_offsets.last().cloned().unwrap();
    let g_pol = fft::fft(&g_hat,
                         rt.pow(last_folding_param as u128),
                         last_eval_offset.pow(last_folding_param as u128));
    let final_deg = params.maxdeg_plus_1 as usize / params.folding_params.iter().product::<usize>();
    if !is_fake {
        assert!(g_pol[(2 * final_deg + 1)..].iter().all(|&x| x == M31::zero()));
    }
    proof.legacy_proof.extend(g_pol[..(2 * final_deg + 1)].iter().flat_map(|c| to_32_be_bytes(c.0 as u128)));
    proof.last_g_pol = g_pol;

    let t_vals = get_pseudorandom_indices(&proof.legacy_proof, 2 * folded_len as u32, params.repetition_params.last().cloned().unwrap(), 0, &mut vec![]);
    let t_shifts = t_vals.iter().map(|&t| t / 2).collect_vec();
    let t_conj = t_vals.iter().map(|&t| t % 2).collect_vec();

    let oracle_branches =
        make_oracle_branches(last_folding_param,
                             params.eval_sizes.last().cloned().unwrap(),
                             &t_shifts, &t_conj, folded_len, &m);
    proof.legacy_proof.extend(oracle_branches.iter().cloned().flatten().flatten().collect_vec());
    proof.last_oracle_branches = oracle_branches;

    proof
}

fn make_oracle_branches(
    folding_param: usize,
    eval_size: usize,
    t_shifts: &[usize],
    t_conj: &[usize],
    folded_len: usize,
    m: &[Hash],
) -> Vec<Vec<Hash>> {
    let mut result = vec![];
    for (&t, &k) in t_shifts.iter().zip(t_conj.iter()) {
        for j in 0..folding_param {
            let branch = mk_branch(&m, t + j * folded_len + k * eval_size);
            result.push(branch)
        }
    }
    result
}

fn verify_low_degree_proof<F: StirField<M31>>(proof: &StirProof, params: &Parameters<F>) -> bool {
    macro_rules! reject_unless_eq {
        ($lhs:expr,$rhs:expr) => { if $lhs != $rhs { return false; } };
    }
    fn to_i128(vs: &[M31]) -> Vec<i128> {
        vs.iter().map(|v| v.0 as i128).collect_vec()
    }

    reject_unless_eq!(params.folding_params.len(), params.eval_offsets.len());
    reject_unless_eq!(params.folding_params.len(), params.eval_sizes.len());
    reject_unless_eq!(params.folding_params.len(), params.repetition_params.len());

    let mut rt = params.root_of_unity;

    reject_unless_eq!(rt.pow(params.eval_sizes[0] as u128), F::one());
    reject_unless_eq!(rt.pow(params.eval_sizes[0] as u128 / 2), F::one() * -F::one());

    let mut proof_pos = 0;
    let mut m_root = &proof.merkle_root;
    proof_pos += 32;
    let mut r_fold =
        params.prim_root.pow(
            (BigUint::from_bytes_be(m_root) % BigUint::from(MODULUS + 1)).to_u128().unwrap(),
        );

    let mut pol: Option<(Vec<M31>, Vec<M31>)> = None;
    let mut rs: Option<Vec<F>> = None;
    let mut zpol: Option<(Vec<M31>, Vec<M31>)> = None;
    let mut r_comb: Option<F> = None;

    for i in 1..params.folding_params.len() {
        reject_unless_eq!(params.eval_sizes[i - 1] % params.folding_params[i - 1], 0);
        let folded_len = params.eval_sizes[i - 1] / params.folding_params[i - 1];
        reject_unless_eq!(params.eval_sizes[i] % folded_len, 0);
        let expand_factor = params.eval_sizes[i] / folded_len;
        reject_unless_eq!(params.eval_sizes[i - 1] % params.eval_sizes[i], 0);
        let eval_size_scale = params.eval_sizes[i - 1] / params.eval_sizes[i];

        let rt_new = rt.pow(eval_size_scale as u128);
        let rt2 = rt_new.pow(expand_factor as u128);
        let p_offset = params.eval_offsets[i - 1].pow(params.folding_params[i - 1] as u128);

        let layer = &proof.layers[i - 1];
        proof_pos += 32;

        let r_outs = get_pseudorandom_element_outside_coset_circle(
            &proof.legacy_proof[..proof_pos], params.prim_root,
            params.eval_sizes[i], params.eval_offsets[i], params.ood_rep,
        );
        let betas: &[M31] = &layer.betas;
        proof_pos += params.ood_rep * 32;

        let r_fold_new = params.prim_root.pow(get_pseudorandom_indices(&proof.legacy_proof[..proof_pos], MODULUS + 1, 1, 0, &mut vec![])[0] as u128);
        let r_comb_new = params.prim_root.pow(get_pseudorandom_indices(&proof.legacy_proof[..proof_pos], MODULUS + 1, 1, 1, &mut vec![])[0] as u128);

        let t_vals = get_pseudorandom_indices(&proof.legacy_proof[..proof_pos], 2 * folded_len as u32, params.repetition_params[i - 1], 2, &mut vec![]);
        let t_shifts = t_vals.iter().map(|&t| t / 2).collect_vec();
        let t_conj = t_vals.iter().map(|&t| t % 2).collect_vec();
        let rs_new: Vec<F> = r_outs.iter().cloned()
            .chain(t_shifts.iter().zip(t_conj.iter())
                .map(|(&t, &k)| conjugate_with_parity(p_offset * rt2.pow(t as u128), k)))
            .collect_vec();

        // should we change to storing the log since it should always be a power of 2?
        let branch_len = (params.eval_sizes[i - 1] as f64).log2().ceil() as usize + 2;
        let num_branches = params.folding_params[i - 1] * params.repetition_params[i - 1];
        proof_pos += branch_len * num_branches * 32;
        let oracle_branches: &[Vec<Hash>] = &layer.oracle_branches;

        let rt2 = rt.pow(folded_len as u128);
        let xs2s = {
            let mut res = vec![get_power_cycle(rt2, F::one())];
            res.push(res[0].clone());
            (&mut res[1][1..]).reverse();
            res
        };

        let mut g_hat = vec![];
        for k in 0_usize..params.repetition_params[i - 1] {
            let mut vals = vec![];
            let x0 = conjugate_with_parity(rt.pow(t_shifts[k] as u128) * params.eval_offsets[i - 1], t_conj[k]);
            for j in 0_usize..params.folding_params[i - 1] {
                let branch = &oracle_branches[k * params.folding_params[i - 1] + j];
                let ind = t_shifts[k] + j * folded_len;
                reject_unless_eq!(verify_branch(m_root, ind + t_conj[k] * params.eval_sizes[i - 1], branch), true);

                let val: M31 = BigUint::from_bytes_be(&branch[0]).to_u32().unwrap().into();

                let new_val = match pol {
                    Some(ref pol) => {
                        let rs = rs.as_ref().unwrap();
                        let r_comb = r_comb.unwrap();
                        let zpol = zpol.as_ref().unwrap();

                        let x = conjugate_with_parity(rt.pow(ind as u128) * params.eval_offsets[i - 1], t_conj[k]);

                        (val - poly_utils::eval_circ_poly_at(pol, &x)) / poly_utils::eval_circ_poly_at(zpol, &x) * poly_utils::geom_sum((x * r_comb).x(), rs.len() as u64)
                    }
                    None => val,
                };
                vals.push(new_val);
            }
            let new_g_hat = poly_utils::eval_circ_poly_at(&poly_utils::circ_lagrange_interp(&xs2s[t_conj[k]], &vals, true), &(r_fold * x0.complex_conjugate()));
            g_hat.push(new_g_hat);
        }

        rt = rt_new;
        m_root = &layer.merkle_root;
        r_fold = r_fold_new;
        r_comb = Some(r_comb_new);
        zpol = Some(poly_utils::circ_zpoly(&rs_new, None));
        let g_rs = betas.iter().cloned().chain(g_hat.iter().cloned()).collect_vec();
        pol = Some(poly_utils::circ_lagrange_interp(&rs_new, &g_rs, false));
        rs = Some(rs_new);
    }

    let final_deg = params.maxdeg_plus_1 as usize / params.folding_params.iter().product::<usize>();
    let g_pol: Vec<M31> = proof.last_g_pol.clone();
    proof_pos += (2 * final_deg + 1) * 32;

    let last_eval_size = *params.eval_sizes.last().unwrap();
    let last_eval_offset = *params.eval_offsets.last().unwrap();
    let last_folding_param = *params.folding_params.last().unwrap();
    let last_repetition_param = *params.repetition_params.last().unwrap();

    reject_unless_eq!(last_eval_size % last_folding_param, 0);
    let folded_len = last_eval_size / last_folding_param;
    let rt2 = rt.pow(last_folding_param as u128);

    let t_vals = get_pseudorandom_indices(&proof.legacy_proof[..proof_pos], 2 * folded_len as u32, params.repetition_params.last().cloned().unwrap(), 0, &mut vec![]);
    let t_shifts = t_vals.iter().map(|&t| t / 2).collect_vec();
    let t_conj = t_vals.iter().map(|&t| t % 2).collect_vec();

    // Warning: code below copy-and-pasted from inside main loop
    let oracle_branches = &proof.last_oracle_branches;

    let rs = rs.as_ref().unwrap();

    let rt3 = rt.pow(folded_len as u128);
    let xs2s = {
        let mut res = vec![get_power_cycle(rt3, F::one())];
        res.push(res[0].clone());
        (&mut res[1][1..]).reverse();
        res
    };

    for k in 0_usize..last_repetition_param {
        let mut vals = vec![];
        let x0 = conjugate_with_parity(rt.pow(t_shifts[k] as u128) * last_eval_offset, t_conj[k]);
        for j in 0..last_folding_param {
            let branch = &oracle_branches[k * last_folding_param + j];
            let ind = t_shifts[k] + j * folded_len;
            reject_unless_eq!(verify_branch(m_root, ind + t_conj[k] * last_eval_size, branch), true);
            let val: M31 = BigUint::from_bytes_be(&branch[0]).to_u32().unwrap().into();

            let x = conjugate_with_parity(rt.pow(ind as u128) * last_eval_offset, t_conj[k]);
            vals.push((
                (val - poly_utils::eval_circ_poly_at(&pol.as_ref().unwrap(), &x)) / poly_utils::eval_circ_poly_at(zpol.as_ref().unwrap(), &x)
            ) * poly_utils::geom_sum((x * r_comb.unwrap()).x(), rs.len() as u64));
        }

        reject_unless_eq!(poly_utils::eval_circ_poly_at(&poly_utils::circ_lagrange_interp(&xs2s[t_conj[k]], &vals, true), &(r_fold * x0.complex_conjugate())),
                          fft::fft_inv(&g_pol, F::one(),
                                       conjugate_with_parity(rt2.pow(t_shifts[k] as u128) * last_eval_offset.pow(last_folding_param as u128),
                                                             t_conj[k]))[0]);
    }

    println!("STIR proof verified");
    true
}

fn get_pseudorandom_element_outside_coset_circle<B, F: StirField<B>>(
    seed: &[u8],
    prim_root: F,
    coset_size: usize,
    shift: F,
    count: usize,
) -> Vec<F> {
    let adjust = 1;
    let m2 = MODULUS + 1;
    assert_eq!((MODULUS as i128 + adjust) % (coset_size as i128), 0);
    let cofactor = (MODULUS as i128 + adjust) / (coset_size as i128);
    let mut exclude = vec![];
    let mut ans = vec![];
    let mut start = 0;
    while ans.len() < count {
        let r = get_pseudorandom_indices(seed, m2 - coset_size as u32, count, start, &mut exclude)[0];
        start += 1;
        exclude.push(r);
        let t = r + 1 + r / (cofactor as usize - 1);

        let val = prim_root.pow(t as u128) * shift;
        if (val * shift).pow(coset_size as u128) != F::one() {
            ans.push(val);
        }
    }
    ans
}

fn conjugate_with_parity<F: ComplexConjugate>(f: F, parity: usize) -> F {
    if parity % 2 == 0 {
        f
    } else {
        f.complex_conjugate()
    }
}

#[cfg(test)]
mod tests {
    use crate::cm31;
    use crate::core::fields::FieldExpOps;
    use super::super::fft::fft_inv;
    use super::*;

    #[test]
    fn test_circle_stir() {
        let prim_root = cm31!(311014874, 1584694829);
        let root_of_unity = prim_root.pow(((MODULUS + 1) / 2_u32.pow(10 + 2)) as u128);
        let log_d = 10;
        let params = generate_parameters::<M31, CM31>(prim_root, root_of_unity, prim_root, log_d, 128, 4, 1.0, 3);
        println!("{:?}", params);

        // Pure STIR tests
        let poly: Vec<M31> = (0..2_u32.pow(log_d + 1)).map(|v| v.into()).collect();
        let evaluations = fft_inv(&poly, root_of_unity, prim_root);
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

    fn generate_parameters<B, F: StirField<B>>(
        prim_root: F,
        root_of_unity: F,
        init_offset: F,
        log_d: u32,
        _security_param: u64,
        log_stopping_degree: u64,
        proximity_param: f64,
        log_folding_param: u32,
    ) -> Parameters<F> {
        let m = ((log_d as u64 - log_stopping_degree) / log_folding_param as u64) as usize;
        let size_l = get_power_cycle(root_of_unity, F::one()).len();
        assert!(size_l.is_power_of_two());
        assert!(proximity_param > 0.0 && proximity_param <= 1.0);

        let mut eval_sizes = vec![size_l];
        let mut eval_offsets = vec![init_offset];
        let mut rt = root_of_unity;

        for _ in 0..m {
            eval_sizes.push(eval_sizes.last().unwrap() / 2);
            eval_offsets.push(rt * init_offset);
            rt = rt.square();
        }

        Parameters {
            root_of_unity,
            maxdeg_plus_1: 2_u32.pow(log_d),
            folding_params: vec![2_usize.pow(log_folding_param); m + 1],
            eval_offsets,
            eval_sizes,
            repetition_params: vec![3; m + 1],
            ood_rep: 1,
            prim_root,
        }
    }
}
