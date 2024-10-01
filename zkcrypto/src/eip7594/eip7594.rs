use crate::kzg_proofs::KZGSettings;
use crate::kzg_types::{Fr, G1, C_KZG_RET};
use zkcrypto::constants;
use zkcrypto::bytes::bytes_to_bls_field;
use kzg::eip_4844::BYTES_PER_FIELD_ELEMENT;
use ckzg::finite::{blst_fr_mul, blst_fr_add};

// remaining helper functions
// bit_reversal_permutation
// fr_ifft
// get_inv_coset_shift_for_cell
// shift_poly
// g1_lincomb_fast
fn compute_commitment_to_aggregated_interpolation_poly(
    commitment_out: &mut G1,
    r_powers: &[Fr],
    cell_indices: &[u64],
    cells: &[Cell],
    num_cells: u64,
    settings: &KZGSettings
) -> C_KZG_RET {
    let mut ret: C_KZG_RET;

    ////////////////////////////////////////////////////////////////////////////////////////////////
    // Array allocations
    ////////////////////////////////////////////////////////////////////////////////////////////////

    let mut is_cell_used = vec![false; CELLS_PER_EXT_BLOB];
    let mut aggregated_column_cells = vec![Fr::zero(); FIELD_ELEMENTS_PER_EXT_BLOB];
    let mut column_interpolation_poly = vec![Fr::zero(); FIELD_ELEMENTS_PER_CELL];
    let mut aggregated_interpolation_poly = vec![Fr::zero(); FIELD_ELEMENTS_PER_CELL];

    ////////////////////////////////////////////////////////////////////////////////////////////////
    // Aggregate cells from the same column
    ////////////////////////////////////////////////////////////////////////////////////////////////

    // Zeroed out columns
    for i in 0..CELLS_PER_EXT_BLOB {
        for j in 0..FIELD_ELEMENTS_PER_CELL {
            let index = i * FIELD_ELEMENTS_PER_CELL + j;
            aggregated_column_cells[index] = Fr::zero();
        }
    }

    // Vertically collapse cells
    for cell_index in 0..num_cells as usize {
        let column_index = cell_indices[cell_index] as usize;

        for fr_index in 0..FIELD_ELEMENTS_PER_CELL {
            let offset = fr_index * BYTES_PER_FIELD_ELEMENT;
            let original_fr = bytes_to_bls_field(&cells[cell_index].bytes[offset..offset + BYTES_PER_FIELD_ELEMENT])?;
            let scaled_fr = blst_fr_mul(&original_fr, &r_powers[cell_index]);

            let array_index = column_index * FIELD_ELEMENTS_PER_CELL + fr_index;
            aggregated_column_cells[array_index] = blst_fr_add(&aggregated_column_cells[array_index], &scaled_fr);
        }
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////
    // Determine which cells are used
    ////////////////////////////////////////////////////////////////////////////////////////////////

    for i in 0..num_cells as usize {
        is_cell_used[cell_indices[i] as usize] = true;
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////
    // Compute interpolation polynomials using the aggregated cells
    ////////////////////////////////////////////////////////////////////////////////////////////////

    // Zero out the polynomial
    for i in 0..FIELD_ELEMENTS_PER_CELL {
        aggregated_interpolation_poly[i] = Fr::zero();
    }

    // Interpolate each column
    for i in 0..CELLS_PER_EXT_BLOB {
        if !is_cell_used[i] {
            continue;
        }

        let index = i * FIELD_ELEMENTS_PER_CELL;
        ret = bit_reversal_permutation(&mut aggregated_column_cells[index..index + FIELD_ELEMENTS_PER_CELL])?;

        ret = fr_ifft(
            &mut column_interpolation_poly,
            &aggregated_column_cells[index..index + FIELD_ELEMENTS_PER_CELL],
            FIELD_ELEMENTS_PER_CELL,
            settings,
        )?;

        let inv_coset_factor = get_inv_coset_shift_for_cell(i, settings);
        shift_poly(&mut column_interpolation_poly, FIELD_ELEMENTS_PER_CELL, &inv_coset_factor);

        for k in 0..FIELD_ELEMENTS_PER_CELL {
            aggregated_interpolation_poly[k] = blst_fr_add(&aggregated_interpolation_poly[k], &column_interpolation_poly[k]);
        }
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////
    // Commit to the aggregated interpolation polynomial
    ////////////////////////////////////////////////////////////////////////////////////////////////

    ret = g1_lincomb_fast(
        commitment_out,
        &settings.g1_values_monomial,
        &aggregated_interpolation_poly,
        FIELD_ELEMENTS_PER_CELL,
    )?;

    Ok(ret)
}
