#[cfg(test)]
mod tests {
    use kzg_bench::tests::eip_7594::{
        test_vectors_verify_cell_kzg_proof_batch,
    };
    use rust_kzg_zkcrypto::{
        eip_4844::load_trusted_setup_filename_rust,
        eip_7594::{
            verify_cell_kzg_proof_batch,
        },
        kzg_proofs::{
            KZGSettings,
            FFTSettings
        },
        kzg_types::{
            ZFp,
            ZFr,
            ZG1,
            ZG1Affine,
            ZG2,
        },
        poly::PolyData,
    };

    #[test]
    pub fn test_vectors_verify_cell_kzg_proof_batch_() {
        test_vectors_verify_cell_kzg_proof_batch::<
            ZFr,
            ZG1,
            ZG2,
            PolyData,
            FFTSettings,
            KZGSettings,
            ZFp,
            ZG1Affine,
        >(
            &load_trusted_setup_filename_rust,
            &verify_cell_kzg_proof_batch,
        );
    }
}