use kzg:: { eip_4844::{ compute_powers, FIELD_ELEMENTS_PER_CELL, FIELD_ELEMENTS_PER_EXT_BLOB, BYTES_PER_COMMITMENT, CELLS_PER_EXT_BLOB, RANDOM_CHALLENGE_KZG_CELL_BATCH_DOMAIN, BYTES_PER_CELL, BYTES_PER_PROOF, hash, hash_to_bls_field, FIELD_ELEMENTS_PER_BLOB }, Fr, G1, G2, common_utils::{reverse_bit_order, reverse_bits_limited}, G1Mul };
use crate::{kzg_proofs::{KZGSettings, pairings_verify, FFTSettings}, kzg_types::{ZFr, ZG1, ZG2}, fft_g1::{g1_linear_combination, fft_g1_fast}, fft::fft_fr_fast};

static CELL_INDICES_RBL: [usize; CELLS_PER_EXT_BLOB] = [
    0x00, 0x40, 0x20, 0x60, 0x10, 0x50, 0x30, 0x70, 0x08, 0x48, 0x28, 0x68, 0x18, 0x58, 0x38, 0x78,
    0x04, 0x44, 0x24, 0x64, 0x14, 0x54, 0x34, 0x74, 0x0c, 0x4c, 0x2c, 0x6c, 0x1c, 0x5c, 0x3c, 0x7c,
    0x02, 0x42, 0x22, 0x62, 0x12, 0x52, 0x32, 0x72, 0x0a, 0x4a, 0x2a, 0x6a, 0x1a, 0x5a, 0x3a, 0x7a,
    0x06, 0x46, 0x26, 0x66, 0x16, 0x56, 0x36, 0x76, 0x0e, 0x4e, 0x2e, 0x6e, 0x1e, 0x5e, 0x3e, 0x7e,
    0x01, 0x41, 0x21, 0x61, 0x11, 0x51, 0x31, 0x71, 0x09, 0x49, 0x29, 0x69, 0x19, 0x59, 0x39, 0x79,
    0x05, 0x45, 0x25, 0x65, 0x15, 0x55, 0x35, 0x75, 0x0d, 0x4d, 0x2d, 0x6d, 0x1d, 0x5d, 0x3d, 0x7d,
    0x03, 0x43, 0x23, 0x63, 0x13, 0x53, 0x33, 0x73, 0x0b, 0x4b, 0x2b, 0x6b, 0x1b, 0x5b, 0x3b, 0x7b,
    0x07, 0x47, 0x27, 0x67, 0x17, 0x57, 0x37, 0x77, 0x0f, 0x4f, 0x2f, 0x6f, 0x1f, 0x5f, 0x3f, 0x7f,
];

pub fn recover_cells_and_kzg_proofs(
    recovered_cells: &mut [[ZFr; FIELD_ELEMENTS_PER_CELL]],
    recovered_proofs: Option<&mut [ZG1]>,
    cell_indices: &[usize],
    cells: &[[ZFr; FIELD_ELEMENTS_PER_CELL]],
    s: &KZGSettings
) -> Result<(), String> {
    if cells.len() > CELLS_PER_EXT_BLOB {
        return Err(String::from("Cell length cannot be larger than CELLS_PER_EXT_BLOB"))
    }

    if cells.len() != cell_indices.len() {
        return Err(
            "Cell indicies mismatch - cells length must be equal to cell indicies length"
                .to_string(),
        );
    }

    if cells.len() < CELLS_PER_EXT_BLOB / 2 {
        return Err(String::from("Impossible to recover - cells length cannot be less than CELLS_PER_EXT_BLOB / 2"));
    }

    for cell_index in cell_indices {
        if *cell_index >= CELLS_PER_EXT_BLOB {
            return Err(String::from("Cell index cannot be larger than CELLS_PER_EXT_BLOB"));
        }
    }

    for cell in recovered_cells.iter_mut() {
        for fr in cell {
            *fr = ZFr::null();
        }
    }
    // ERROR
    for i in 0..cells.len() {
        let index = cell_indices[i]; // * FIELD_ELEMENTS_PER_CELL;
        for j in 0..FIELD_ELEMENTS_PER_CELL {
            if !recovered_cells[index][j].is_null() {
                return Err(String::from("Invalid output cell"));
            }
        }

        recovered_cells[index] = cells[i];
    }

    if cells.len() != CELLS_PER_EXT_BLOB {
        recover_cells(recovered_cells.as_flattened_mut(), cell_indices, s)?;
    }
    
    if let Some(recovered_proofs) = recovered_proofs {
        let mut poly = vec![ZFr::default(); FIELD_ELEMENTS_PER_EXT_BLOB];

        poly_lagrange_to_monomial(&mut poly, recovered_cells.as_flattened(), &s.fs)?;

        compute_fk20_proofs(recovered_proofs, &poly, FIELD_ELEMENTS_PER_BLOB, s)?;
    
        reverse_bit_order(recovered_proofs)?;
    }

    return Ok(());
}

fn compute_fk20_proofs(
    proofs: &mut [ZG1],
    poly: &[ZFr],
    n: usize,
    s: &KZGSettings,
) -> Result<(), String> {
    let k = n / FIELD_ELEMENTS_PER_CELL;
    let k2 = k * 2;

    let mut coeffs = vec![vec![ZFr::default(); k]; k2];
    let mut h_ext_fft = vec![ZG1::identity(); k2];
    let mut toeplitz_coeffs = vec![ZFr::default(); k2];
    let mut toeplitz_coeffs_fft = vec![ZFr::default(); k2];

    for i in 0..FIELD_ELEMENTS_PER_CELL {
        toeplitz_coeffs_stride(&mut toeplitz_coeffs, poly, n, i, FIELD_ELEMENTS_PER_CELL)?;
        fr_fft(&mut toeplitz_coeffs_fft, &toeplitz_coeffs, &s.fs)?;
        for j in 0..k2 {
            coeffs[j][i] = toeplitz_coeffs_fft[j];
        }
    }

    for i in 0..k2 {
        g1_linear_combination(
            &mut h_ext_fft[i],
            &s.x_ext_fft_columns[i],
            &coeffs[i],
            FIELD_ELEMENTS_PER_CELL,
            None,
        );
    }

    let mut h = vec![ZG1::identity(); k2];
    g1_ifft(&mut h, &h_ext_fft, s, &k2)?;

    for i in k..k2 {
        h[i] = ZG1::identity();
    }

    g1_fft(proofs, &h, s, &k2)?;

    Ok(())
}

fn g1_fft(
    out: &mut [ZG1], 
    input: &[ZG1], 
    s: &KZGSettings, 
    k2: &usize
) -> Result<(), String> {
    if input.len() > FIELD_ELEMENTS_PER_EXT_BLOB || !input.len().is_power_of_two() {
        return Err("Invalid input length".to_string());
    }

    let roots_stride = FIELD_ELEMENTS_PER_EXT_BLOB / input.len();
    fft_g1_fast(out, input, 1, &s.fs.roots_of_unity, roots_stride, *k2);

    Ok(())
}

fn g1_ifft(
    out: &mut [ZG1], 
    input: &[ZG1],
    s: &KZGSettings,
    k2: &usize,
) -> Result<(), String> {
    if input.len() > FIELD_ELEMENTS_PER_EXT_BLOB || !input.len().is_power_of_two() {
        return Err("Invalid input length".to_string());
    }

    let stride = FIELD_ELEMENTS_PER_EXT_BLOB / input.len();
    fft_g1_fast(out, &input, 1, &s.fs.reverse_roots_of_unity, stride, *k2);

    let inv_len = ZFr::from_u64(input.len() as u64).eucl_inverse();
    for i in 0..input.len() {
        out[i] = out[i].mul(&inv_len);
    }

    Ok(())
}

fn toeplitz_coeffs_stride(
    out: &mut [ZFr],
    input: &[ZFr],
    n: usize,
    offset: usize,
    stride: usize,
) -> Result<(), String> {
    if stride == 0 {
        return Err(String::from("Stride cannot be zero"));
    }

    let k = n / stride;
    let k2 = k * 2;

    out[0] = input[n - 1 - offset];
    {
        let mut i = 1;
        while i <= k + 1 && i < k2 {
            out[i] = ZFr::zero();
            i += 1;
        }
    };

    {
        let mut i = k + 2;
        let mut j = 2 * stride - offset - 1;
        while i < k2 {
            out[i] = input[j];
            i += 1;
            j += stride;
        }
    };

    Ok(())
}

fn poly_lagrange_to_monomial(
    output: &mut [ZFr],
    lagrange_poly: &[ZFr],
    s: &FFTSettings,
) -> Result<(), String> {
    let mut poly = lagrange_poly.to_vec();

    reverse_bit_order(&mut poly)?;

    fr_ifft(output, &poly, s)?;

    Ok(())
}

fn recover_cells(
    output: &mut [ZFr],
    cell_indices: &[usize],
    s: &KZGSettings,
) -> Result<(), String> {
    let mut missing_cell_indices: Vec<usize> = Vec::new();

    let mut cells_brp = output.to_vec();
    reverse_bit_order(&mut cells_brp)?;

    for i in 0..CELLS_PER_EXT_BLOB {
        if !cell_indices.contains(&i) {
            missing_cell_indices.push(reverse_bits_limited(CELLS_PER_EXT_BLOB, i));
        }
    }

    if missing_cell_indices.len() > CELLS_PER_EXT_BLOB / 2 {
        return Err(String::from("Not enough cells"));
    }

    let missing_cell_indices = &missing_cell_indices[..];
    let vanishing_poly_coeff = vanishing_polynomial_for_missing_cells(missing_cell_indices, s)?;

    let mut vanishing_poly_eval = vec![ZFr::default(); FIELD_ELEMENTS_PER_EXT_BLOB];

    fr_fft(&mut vanishing_poly_eval, &vanishing_poly_coeff, &s.fs)?;

    let mut extended_evaluation_times_zero: Vec<ZFr> = Vec::with_capacity(FIELD_ELEMENTS_PER_EXT_BLOB);

    for i in 0..FIELD_ELEMENTS_PER_EXT_BLOB {
        if cells_brp[i].is_null() {
            extended_evaluation_times_zero.push(ZFr::zero());
        } else {
            extended_evaluation_times_zero.push(cells_brp[i].mul(&vanishing_poly_eval[i]));
        }
    }

    let mut extended_evaluation_times_zero_coeffs = vec![ZFr::default(); FIELD_ELEMENTS_PER_EXT_BLOB];

    fr_ifft(&mut extended_evaluation_times_zero_coeffs, &extended_evaluation_times_zero, &s.fs)?;

    let mut extended_evaluations_over_coset = coset_fft(&extended_evaluation_times_zero_coeffs, &s)?;

    let vanishing_poly_over_coset = coset_fft(&vanishing_poly_coeff, &s)?;

    for i in 0..FIELD_ELEMENTS_PER_EXT_BLOB {
        extended_evaluations_over_coset[i] = extended_evaluations_over_coset[i].div(&vanishing_poly_over_coset[i])?;
    }

    let reconstructed_poly_coeff = coset_ifft(&extended_evaluations_over_coset, s)?;

    fr_fft(output, &reconstructed_poly_coeff, &s.fs)?;

    reverse_bit_order(output)?;

    return Ok(())
}

fn coset_ifft(input: &[ZFr], s: &KZGSettings) -> Result<Vec<ZFr>, String> {
    if input.len() == 0 {
        return Err("Invalid input length".to_string());
    }

    let mut output = vec![ZFr::default(); input.len()];
    fr_ifft(&mut output, input, &s.fs)?;

    shift_poly(&mut output, &ZFr::one().div(&ZFr::from_u64(7))?);

    Ok(output)
}

fn coset_fft(input: &[ZFr], s: &KZGSettings) -> Result<Vec<ZFr>, String> {
    if input.len() == 0 {
        return Err("Invalid input length".to_string());
    }

    let mut in_shifted = input.to_vec();
    // TODO: move 7 to constant
    shift_poly(&mut in_shifted, &ZFr::from_u64(7));

    let mut output = vec![ZFr::default(); input.len()];
    fr_fft(&mut output, &in_shifted, &s.fs)?;

    Ok(output)
}

fn fr_ifft(output: &mut [ZFr], input: &[ZFr], s: &FFTSettings) -> Result<(), String> {
    let stride = FIELD_ELEMENTS_PER_EXT_BLOB / input.len();

    fft_fr_fast(output, &input, 1, &s.reverse_roots_of_unity, stride);

    let inv_len = ZFr::from_u64(input.len().try_into().unwrap()).inverse();
    for el in output {
        *el = el.mul(&inv_len);
    }

    Ok(())
}


fn fr_fft(output: &mut [ZFr], input: &[ZFr], s: &FFTSettings) -> Result<(), String> {
    let roots_stride = FIELD_ELEMENTS_PER_EXT_BLOB / input.len();
    fft_fr_fast(output, &input, 1, &s.roots_of_unity, roots_stride);

    Ok(())
}


pub fn vanishing_polynomial_for_missing_cells(
    missing_cell_indices: &[usize],
    s: &KZGSettings,
) -> Result<Vec<ZFr>, String> {
    if missing_cell_indices.len() == 0 || missing_cell_indices.len() >= CELLS_PER_EXT_BLOB {
        return Err(String::from("Invalid missing cell indices count"));
    }

    const STRIDE: usize = FIELD_ELEMENTS_PER_EXT_BLOB / CELLS_PER_EXT_BLOB;

    // let roots = missing_cell_indices
    //     .iter()
    //     .map(|i| s.fs.roots_of_unity[*i * STRIDE])
    //     .collect::<Vec<_>>();

    let mut roots = Vec::new();
    for i in missing_cell_indices {
        let a = *i;
        let b = *i * STRIDE;
        roots.push(s.fs.roots_of_unity[*i * STRIDE]);
    }

    let short_vanishing_poly = compute_vanishing_polynomial_from_roots(&roots)?;

    let mut vanishing_poly = vec![ZFr::zero(); FIELD_ELEMENTS_PER_EXT_BLOB];
    
    for i in 0..short_vanishing_poly.len() {
        vanishing_poly[i * FIELD_ELEMENTS_PER_CELL] = short_vanishing_poly[i];
    }

    Ok(vanishing_poly)
}   

fn compute_vanishing_polynomial_from_roots(
    roots: &[ZFr],
) -> Result<Vec<ZFr>, String> {
    if roots.len() == 0 {
        return Err(String::from("Roots cannot be empty"));
    }

    let mut poly = Vec::new();
    poly.push(roots[0].negate());

    for i in 1..roots.len() {
        let neg_root = roots[i].negate();

        poly.push(neg_root.add(&poly[i - 1]));

        for j in (1..i).rev() {
            poly[j] = poly[j].mul(&neg_root).add(&poly[j - 1]);
        }
        poly[0] = poly[0].mul(&neg_root);
    }

    poly.push(ZFr::one());

    Ok(poly)    
}

pub fn verify_cell_kzg_proof_batch(
    commitments_bytes: &[ZG1],
    cell_indices: &[usize],
    cells: &[[ZFr; FIELD_ELEMENTS_PER_CELL]],
    proofs_bytes: &[ZG1],
    s: &KZGSettings,
) -> Result<bool, String> {
    /* Exit early if we are given zero cells */
    if cells.len() == 0 {
        return Ok(true);
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////
    // Sanity checks
    ////////////////////////////////////////////////////////////////////////////////////////////////
    for cell_index in cell_indices {
        /* Make sure column index is valid */
        if *cell_index >= CELLS_PER_EXT_BLOB {
            return Err(String::from("The supplied data is invalid"));
        }
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////
    // Deduplicate commitments
    ////////////////////////////////////////////////////////////////////////////////////////////////
    // let mut unique_commitments = vec![ZFp(Fp::default()); *num_cells];
    let mut unique_commitments = commitments_bytes.to_vec();
    let mut commitment_indices = vec![0usize; cells.len()];
    let mut num_commitments = commitments_bytes.len();
    
    unique_commitments.extend_from_slice(commitments_bytes);
    deduplicate_commitments(&mut unique_commitments, &mut commitment_indices, &mut num_commitments);
    
    // let mut r_powers = vec![ZFr::zero(); *num_cells];
    // let mut proofs_g1 = vec![ZG1::zero(); *num_cells];
    let unique_commitments = &unique_commitments[0..num_commitments];
    let r_powers = compute_r_powers_for_verify_cell_kzg_proof_batch(
        &unique_commitments,
        commitment_indices.len(),
        &commitment_indices,
        cell_indices,
        cells,
        proofs_bytes,
        cells.len(),
    )?;

    let mut proof_lincomb = ZG1::default();
    g1_linear_combination(&mut proof_lincomb, &proofs_bytes, &r_powers, cells.len(), None);

    let final_g1_sum = compute_weighted_sum_of_commitments(
        &unique_commitments, 
        &commitment_indices, 
        &r_powers,
    );

    let interpolation_poly_commit = compute_commitment_to_aggregated_interpolation_poly(
        &r_powers, 
        &cell_indices, 
        &cells, 
        s,
    );

    match interpolation_poly_commit {
        Ok(res) => final_g1_sum.sub(&res),
        Err(msg) => return Err(msg),
    };


    let weighted_sum_of_proofs = computed_weighted_sum_of_proofs(
        proofs_bytes, 
        &r_powers, 
        cell_indices, 
        &proofs_bytes.len(), 
        s
    )?;

    let final_g1_sum = final_g1_sum.add(&weighted_sum_of_proofs);

    let power_of_s = &s.g2_values_monomial[FIELD_ELEMENTS_PER_CELL];
    
    Ok(pairings_verify(
        &final_g1_sum,
        &ZG2::generator(),
        &proof_lincomb,
        &power_of_s,
    ))
}

// fn fr_ifft(output: &mut [ZFr], input: &[ZFr], n: usize, s: &FFTSettings) -> Result<(), String> {
//     if n == 0 {
//         return Ok(());
//     }

//     if n > FIELD_ELEMENTS_PER_EXT_BLOB || !is_power_of_2(n) {
//         return Err(String::from("Bad arguments"));
//     }

//     let stride = FIELD_ELEMENTS_PER_EXT_BLOB / n;

//     fft_fr_fast(output, &input, 1, &s.reverse_roots_of_unity, stride);

//     let inv_len = ZFr::from_u64(input.len().try_into().unwrap()).inverse();
//     for el in output {
//         *el = el.mul(&inv_len);
//     }

//     Ok(())
// }

fn compute_commitment_to_aggregated_interpolation_poly(
    r_powers: &[ZFr],
    cell_indices: &[usize],
    cells: &[[ZFr; FIELD_ELEMENTS_PER_CELL]],
    s: &KZGSettings,
) -> Result<ZG1, String> {
    let mut aggregated_column_cells =
        vec![ZFr::zero(); CELLS_PER_EXT_BLOB * FIELD_ELEMENTS_PER_CELL];

    for (cell_index, column_index) in cell_indices.iter().enumerate() {
        for fr_index in 0..FIELD_ELEMENTS_PER_CELL {
            let original_fr = cells[cell_index][fr_index];

            let scaled_fr = original_fr.mul(&r_powers[cell_index]);

            let array_index = column_index * FIELD_ELEMENTS_PER_CELL + fr_index;
            aggregated_column_cells[array_index] =
                aggregated_column_cells[array_index].add(&scaled_fr);
        }
    }

    let mut is_cell_used = vec![false; CELLS_PER_EXT_BLOB];

    for cell_index in cell_indices {
        is_cell_used[*cell_index] = true;
    }

    let mut aggregated_interpolation_poly = vec![ZFr::zero(); FIELD_ELEMENTS_PER_CELL];
    let mut column_interpolation_poly = vec![ZFr::default(); FIELD_ELEMENTS_PER_CELL];
    for i in 0..CELLS_PER_EXT_BLOB {
        if !is_cell_used[i] {
            continue;
        }

        let index = i * FIELD_ELEMENTS_PER_CELL;

        reverse_bit_order(&mut aggregated_column_cells[index..(index + FIELD_ELEMENTS_PER_CELL)])?;

        fr_ifft(
            &mut column_interpolation_poly,
            &aggregated_column_cells[index..(index + FIELD_ELEMENTS_PER_CELL)],
            &s.fs,
        )?;

        let inv_coset_factor = get_inv_coset_shift_for_cell(i, s)?;

        shift_poly(&mut column_interpolation_poly, &inv_coset_factor);

        for k in 0..FIELD_ELEMENTS_PER_CELL {
            aggregated_interpolation_poly[k] =
                aggregated_interpolation_poly[k].add(&column_interpolation_poly[k]);
        }
    }

    let mut commitment_out = ZG1::default();
    g1_linear_combination(
        &mut commitment_out,
        &s.g1_values_monomial,
        &aggregated_interpolation_poly,
        FIELD_ELEMENTS_PER_CELL,
        None,
    );

    Ok(commitment_out)
}

fn shift_poly(poly: &mut [ZFr], shift_factor: &ZFr) {
    let mut factor_power = ZFr::one();
    for i in 1..poly.len() {
        factor_power = factor_power.mul(shift_factor);
        poly[i] = poly[i].mul(&factor_power);
    }
}

fn get_inv_coset_shift_for_cell(
    cell_index: usize,
    settings: &KZGSettings,
) -> Result<ZFr, String> {
    /*
     * Get the cell index in reverse-bit order.
     * This index points to this cell's coset factor h_k in the roots_of_unity array.
     */
    let cell_index_rbl = CELL_INDICES_RBL[cell_index];

    if cell_index_rbl > FIELD_ELEMENTS_PER_EXT_BLOB {
        return Err(String::from("Invalid cell index"));
    }
    let inv_coset_factor_idx = FIELD_ELEMENTS_PER_EXT_BLOB - cell_index_rbl;

    /* Get h_k^{-1} using the index */
    if inv_coset_factor_idx >= FIELD_ELEMENTS_PER_EXT_BLOB + 1 {
        return Err(String::from("Invalid cell index"));
    }

    return Ok(settings.fs.roots_of_unity[inv_coset_factor_idx]);
}

fn compute_weighted_sum_of_commitments(
    commitments: &[ZG1],
    commitment_indices: &[usize],
    r_powers: &[ZFr],
) -> ZG1 {
    let mut commitment_weights = vec![ZFr::zero(); commitments.len()];

    for i in 0..r_powers.len() {
        commitment_weights[commitment_indices[i]] =
            commitment_weights[commitment_indices[i]].add(&r_powers[i]);
    }

    let mut sum_of_commitments = ZG1::default();
    g1_linear_combination(
        &mut sum_of_commitments,
        &commitments,
        &commitment_weights,
        commitments.len(),
        None,
    );

    sum_of_commitments
}

fn compute_r_powers_for_verify_cell_kzg_proof_batch(
    commitments_bytes: &[ZG1],
    num_commitments: usize,
    commitment_indices: &[usize],
    cell_indices: &[usize],
    cells: &[[ZFr; FIELD_ELEMENTS_PER_CELL]],
    proofs_bytes: &[ZG1],
    num_cells: usize 
) -> Result<Vec<ZFr>, String> {
    let mut bytes: Vec<u8> = Vec::new(); // Vec::new(); is the same as vec![];
    
    /* Calculate the size of the data we're going to hash */
    let input_size = RANDOM_CHALLENGE_KZG_CELL_BATCH_DOMAIN.len() /* The domain separator */ 
    + std::mem::size_of::<u64>()                               /* FIELD_ELEMENTS_PER_CELL */
    + std::mem::size_of::<u64>()                               /* num_commitments */
    + std::mem::size_of::<u64>()                               /* num_cells */
    + (num_commitments * BYTES_PER_COMMITMENT)                 /* commitment_bytes */
    + (num_cells as usize * std::mem::size_of::<u64>())        /* commitment_indices */ // num_cells as usize because "pub const fn size_of<T>() -> usize" and input_size is usize
    + (num_cells as usize * std::mem::size_of::<u64>())        /* cell_indices */
    + (num_cells as usize * BYTES_PER_CELL)                    /* cells */
    + (num_cells as usize * BYTES_PER_PROOF);                  /* proofs_bytes */

    bytes.resize(input_size, 0); //as an alternative

    let mut bytes = vec![0; input_size];
    bytes[..16].copy_from_slice(&RANDOM_CHALLENGE_KZG_CELL_BATCH_DOMAIN);
    bytes[16..24].copy_from_slice(&(FIELD_ELEMENTS_PER_CELL as u64).to_be_bytes());
    bytes[24..32].copy_from_slice(&(commitments_bytes.len() as u64).to_be_bytes());
    bytes[32..40].copy_from_slice(&(cells.len() as u64).to_be_bytes());

    let mut offset = 40;

    for commitment in commitments_bytes {
        /* Copy commitment */
        bytes[offset..offset + BYTES_PER_COMMITMENT].copy_from_slice(&commitment.to_bytes());
        offset += BYTES_PER_COMMITMENT;
    }
        
    for i in 0..num_cells as usize {
        bytes[offset..(offset + 8)].copy_from_slice(&(commitment_indices[i] as u64).to_be_bytes());
        offset += 8;

        bytes[offset..(offset + 8)].copy_from_slice(&(cell_indices[i] as u64).to_be_bytes());
        offset += 8;

        bytes[offset..(offset + BYTES_PER_CELL)].copy_from_slice(
            &cells[i]
                .iter()
                .flat_map(|fr| fr.to_bytes())
                .collect::<Vec<_>>(),
        );
        offset += BYTES_PER_CELL;

        bytes[offset..(offset + BYTES_PER_PROOF)].copy_from_slice(&(proofs_bytes[i].to_bytes()));
        offset += BYTES_PER_PROOF;
    }
    let bytes = &bytes[..];

    if offset != input_size {
        return Err(String::from("Failed to create challenge - invalid length"));
    }

    let r_bytes_ = hash(&bytes);
    let r = hash_to_bls_field(&r_bytes_);
    let r_powers_out = compute_powers(&r, num_cells as usize);

    Ok(r_powers_out)
}

pub fn bytes_of_uint64(out: &mut [u8], mut n: u64) {
    for byte in out.iter_mut().rev().take(8) {
        *byte = (n & 0xff) as u8;
        n >>= 8;
    }
}

fn computed_weighted_sum_of_proofs(
    proofs_g1: &[ZG1],
    r_powers: &[ZFr],
    cell_indices: &[usize],
    num_cells: &usize,
    s: &KZGSettings,
) -> Result<ZG1, String> {
    let mut weighted_powers_of_r = vec![ZFr::zero(); *num_cells];
    // let mut weighted_powers_of_r = Vec::with_capacity(*num_cells);
    
    for i in 0..(*num_cells) {
        let h_k_pow = get_coset_shift_pow_for_cell(&cell_indices[i], &s);
        // weighted_powers_of_r.push(r_powers[i].mul(&h_k_pow));
        match h_k_pow {
            Ok(result) => weighted_powers_of_r.push(r_powers[i].mul(&result)),
            Err(msg) => return Err(msg),
        }
        // unsafe { blst_fr_mul(&mut weighted_powers_of_r[i].to_blst_fr(), &r_powers[i].to_blst_fr(), &h_k_pow.to_blst_fr()); }
    }
    
    // g1_lincomb_fast(weighted_proof_sum_out, &proofs_g1, &weighted_powers_of_r, &num_cells);
    let mut sum_of_commitments = ZG1::default();
    g1_linear_combination(
        &mut sum_of_commitments,
        &proofs_g1,
        &weighted_powers_of_r,
        *num_cells,
        None,
    );

    Ok(sum_of_commitments)
}

fn get_coset_shift_pow_for_cell(
    cell_index: &usize,
    s: &KZGSettings
) -> Result<ZFr, String> {
    let cell_idx_rbl = CELL_INDICES_RBL[*cell_index];
    let h_k_pow_idx = cell_idx_rbl * FIELD_ELEMENTS_PER_CELL;
    
    if h_k_pow_idx >= FIELD_ELEMENTS_PER_EXT_BLOB + 1 {
        return Err(String::from("Index is invalid"));
    }
    
    // (*s).fs.roots_of_unity[h_k_pow_idx as usize]
    // Ok(s.get_roots_of_unity_at(h_k_pow_idx))
    Ok((*s).fs.roots_of_unity[h_k_pow_idx as usize])
}

fn deduplicate_commitments(
    commitments_out: &mut Vec<ZG1>,
    indices_out: &mut Vec<usize>,
    count_out: &mut usize,    
) {
    // Bail early if there are no commitments
    if *count_out == 0 {
        return;
    }
    
    // The first commitment is always new
    indices_out[0] = 0;
        let mut new_count: usize = 1;
        let mut exist = false;

        // Create list of unique commitments & indices to them
        for i in 0..(*count_out) {
            for j in 0..new_count {
                if commitments_out[i] == commitments_out[j] {
                    // This commitment already exists
                    indices_out[i] = j;
                    exist = true;
                    break;
                }
            }
            if !exist {
                // This is a new commitment
                commitments_out[new_count] = commitments_out[i];
                indices_out[i] = new_count;
                new_count += 1;
            }
        }
}
    
// fn g1_lincomb_fast(
//     out: &mut ZG1,
//     p: &[ZG1],
//     coeffs: &[ZFr],
//     len: &usize,
// ) {
//     let min_length_threshold: usize = 8;
    
//     if *len < min_length_threshold {
    //         g1_lincomb_naive(out, p, coeffs, len);
//     }

//     let mut p_filtered = vec![blst_p1::default(); *len];
//     let mut p_affine = vec![blst_p1_affine::default(); *len];
//     let mut scalars = vec![blst_scalar::default(); *len];
//     let scratch_size: usize;
//     unsafe { scratch_size = blst_p1s_mult_pippenger_scratch_sizeof(*len); }
//     let mut scratch = vec![limb_t::default(); scratch_size];

//     for i in 0..(*len) {
//         unsafe { blst_scalar_from_fr(&mut scalars[i], &coeffs[i].to_blst_fr()); }
//     }

//     let mut new_len: usize = 0;
//     for i in 0..(*len) {
    //         if unsafe { blst_p1_is_inf(&p[i].to_blst_p1()) } {
//             p_filtered[new_len] = p[i].to_blst_p1();
//             scalars[new_len] = scalars[i].clone();
//             new_len += 1;
//         }
//     }

//     if new_len < min_length_threshold {
    //         g1_lincomb_naive(out, p, coeffs, len);
//     }

//     // let p_arg: [&Option<&blst_p1>; 2] = [&Some(&p_filtered[0]), &None];
//     let p_arg = [&p_filtered[0], std::ptr::null()];
//     unsafe { blst_p1s_to_affine(&mut p_affine[0], p_arg.as_ptr(), new_len) };

//     // Sus
//     let scalars_arg = [&scalars[0].b[0], std::ptr::null()];
//     let points_arg = [&p_affine[0], std::ptr::null()];
//     unsafe { blst_p1s_mult_pippenger(&mut out.to_blst_p1(), points_arg.as_ptr(), new_len, scalars_arg.as_ptr(), BITS_PER_FIELD_ELEMENT, scratch.as_mut_ptr()); }
// }

// fn g1_lincomb_naive(
//     out: &mut ZG1,
//     p: &[ZG1],
//     coeffs: &[ZFr],
//     len: &usize,
// ) {
    //     let mut tmp = ZG1::default();
//     *out = ZG1::identity();

//     for i in 0..(*len) {
//         g1_mul(&mut tmp, &p[i], &coeffs[i]);
//         let temp = out.clone();
//         unsafe { blst_p1_add_or_double(&mut out.to_blst_p1(), &temp.to_blst_p1(), &tmp.to_blst_p1()) };
//     }
// }

// fn g1_mul(
//     out: &mut ZG1,
//     a: &ZG1,
//     b: &ZFr
// ) {
//     let mut s: blst_scalar = blst_scalar::default();
//     unsafe { 
//         blst_scalar_from_fr(&mut s, &b.to_blst_fr());
//         blst_p1_mult(&mut out.to_blst_p1(), &a.to_blst_p1(), &s.b[0], BITS_PER_FIELD_ELEMENT);
//     }
// }
// compute_r_powers_for_verify_cell_kzg_proof_batch(
//     r_powers,
//     unique_commitments,
//     num_commitments,
//     commitment_indices,
//     cell_indices,
//     cells,
//     proofs_bytes,
//     num_cells
// );

// bytes_to_kzg_proof(&proofs_g1[i], &proofs_bytes[i]);
// proofs_g1 = ZG1::from_bytes(proofs_bytes).unwrap();