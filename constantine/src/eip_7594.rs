
/**
 * Serialize a 64-bit unsigned integer into bytes.
 *
 * @param[out]  out An 8-byte array to store the serialized integer
 * @param[in]   n   The integer to be serialized
 *
 * @remark The output format is big-endian.
 */
// Written by me
 // fn bytes_of_uint64(out: &mut [u8], mut n: u64) { //maybe &mut[u8;8]?
//     for i in (0..8).rev() {
//         out[i] = (n & 0xFF) as u8;
//         n >>= 8;
//     }
// }
// Written by Artiom
pub fn bytes_of_uint64(out: &mut [u8], mut n: u64) {
    for byte in out.iter_mut().rev().take(8) {
        *byte = (n & 0xff) as u8;
        n >>= 8;
    }
}

/**
 * Compute random linear combination challenge scalars for verify_cell_kzg_proof_batch. In this, we
 * must hash EVERYTHING that the prover can control.
 *
 * @param[out]  r_powers_out        The output challenges, length `num_cells`
 * @param[in]   commitments_bytes   The input commitments, length `num_commitments`
 * @param[in]   num_commitments     The number of commitments
 * @param[in]   commitment_indices  The cell commitment indices, length `num_cells`
 * @param[in]   cell_indices        The cell indices, length `num_cells`
 * @param[in]   cells               The cell, length `num_cells`
 * @param[in]   proofs_bytes        The cell proof, length `num_cells`
 * @param[in]   num_cells           The number of cells
 */
fn compute_r_powers_for_verify_cell_kzg_proof_batch(
    r_powers_out: &mut [fr_t],
    commitments_bytes: &[Bytes48],
    num_commitments: usize,
    commitment_indices: &[u64],
    cell_indices: &[u64],
    cells: &[Cell],
    proofs_bytes: &[Bytes48],
    num_cells: u64 
) -> C_KZG_RET {
    let mut ret: C_KZG_RET;
    let mut bytes: Vec<u8> = Vec::new(); // Vec::new(); is the same as vec![];
    let mut r_bytes: Bytes32;
    let mut r: fr_t; // what is this fr_t????????

    /* Calculate the size of the data we're going to hash */
    let input_size: usize = DOMAIN_STR_LENGTH                         /* The domain separator */ 
    + std::mem::size_of::<u64>()                               /* FIELD_ELEMENTS_PER_CELL */
    + std::mem::size_of::<u64>()                               /* num_commitments */
    + std::mem::size_of::<u64>()                               /* num_cells */
    + (num_commitments * BYTES_PER_COMMITMENT)                 /* commitment_bytes */
    + (num_cells as usize * std::mem::size_of::<u64>())        /* commitment_indices */ // num_cells as usize because "pub const fn size_of<T>() -> usize" and input_size is usize
    + (num_cells as usize * std::mem::size_of::<u64>())        /* cell_indices */
    + (num_cells as usize * BYTES_PER_CELL)                    /* cells */
    + (num_cells as usize * BYTES_PER_PROOF);                  /* proofs_bytes */

    /* Allocate space to copy this data into */ 
    //ret = c_kzg_malloc(&mut bytes, input_size); // if there is no c_kzg_malloc, we will use .resize which is safe, hence no need to check ret!
    bytes.resize(input_size, 0); //as an alternative

    // let mut check_if_not_ok: bool = false; // used instead of goto - will not be needed if bytes.resize is used
    // match ret{
    //     C_KZG_RET::C_KZG_OK => check_if_not_ok = true,
    //     C_KZG_RET::C_KZG_BADARGS => check_if_not_ok = false,
    //     C_KZG_RET::C_KZG_ERROR => check_if_not_ok = false,
    //     C_KZG_RET::C_KZG_MALLOC => check_if_not_ok = false,
    // }

    //if check_if_not_ok {
    /* Pointer tracking `bytes` for writing on top of it */
    let mut offset: u8 = bytes;

    /* Ensure that the domain string is the correct length */
    assert_eq!(RANDOM_CHALLENGE_DOMAIN_VERIFY_CELL_KZG_PROOF_BATCH.len(), DOMAIN_STR_LENGTH);

    /* Copy domain separator */ //CHECK memcpy!
    bytes[offset..offset + DOMAIN_STR_LENGTH]
    .copy_from_slice(&RANDOM_CHALLENGE_DOMAIN_VERIFY_CELL_KZG_PROOF_BATCH); //copy_from_slice is safe alternative of memcpy; copy_nonoverlapping is an unsafe one
    offset += DOMAIN_STR_LENGTH;

    /* Copy field elements per cell */
    bytes_of_uint64(&mut bytes[offset..], FIELD_ELEMENTS_PER_CELL);
    offset += std::mem::size_of::<u64>();

    /* Copy number of commitments */
    bytes_of_uint64(&mut bytes[offset..], num_commitments as u64);
    offset += std::mem::size_of::<u64>();

    /* Copy number of cells */
    bytes_of_uint64(&mut bytes[offset..], num_cells as u64);
    offset += std::mem::size_of::<u64>();

    for commitment in commitments_bytes.iter().take(num_commitments) {
        /* Copy commitment */
        bytes[offset..offset + BYTES_PER_COMMITMENT].copy_from_slice(&commitment.bytes);
        offset += BYTES_PER_COMMITMENT;
    }
    //another option
    // for i in 0..num_commitments as usize { 
    //     /* Copy commitment */
    //     bytes[offset..offset + BYTES_PER_COMMITMENT].copy_from_slice(&commitments_bytes[i].bytes);
    //     other_array[i] = some_computation(i); // Using the index for another purpose
    //     offset += BYTES_PER_COMMITMENT;
    // }
        
    for i in 0..num_cells as usize {
        /* Copy row id */
        bytes_of_uint64(&mut bytes[offset..], commitment_indices[i]);
        offset += std::mem::size_of::<u64>();

        /* Copy column id */
        bytes_of_uint64(&mut bytes[offset..], cell_indices[i]);
        offset += std::mem::size_of::<u64>();

        /* Copy cell */
        bytes[offset..offset + BYTES_PER_CELL].copy_from_slice(&cells[i].bytes);
        offset += BYTES_PER_CELL;

        /* Copy proof */
        bytes[offset..offset + BYTES_PER_PROOF].copy_from_slice(&proofs_bytes[i].bytes);
        offset += BYTES_PER_PROOF;
    }

    /* Now let's create the challenge! */
    blst_sha256(&mut r_bytes.bytes, &bytes, input_size);
    hash_to_bls_field(&mut r, &r_bytes);

    /* Raise power of r for each cell */
    compute_powers(r_powers_out, &r, num_cells as usize);

    /* Make sure we wrote the entire buffer */
    assert_eq!(offset, input_size);
    //}
    ret = C_KZG_RET::C_KZG_OK; //Good, if we don't use c_kzg_malloc
    ret
}
