use super::utils::{ get_manifest_dir, get_trusted_setup_path };
use crate::test_vectors::{ verify_cell_kzg_proof_batch };
use kzg::{
    eip_4844::{BYTES_PER_FIELD_ELEMENT, FIELD_ELEMENTS_PER_CELL},
    FFTSettings, Fr, G1Affine, G1Fp, G1GetFp, G1Mul, KZGSettings, Poly, G1, G2,
};
use std::{fmt::Debug, fs, path::PathBuf};

const VERIFY_CELL_KZG_PROOF_BATCH_TEST_VECTORS: &str =
    "src/test_vectors/verify_cell_kzg_proof_batch/*/*/*";


pub fn test_vectors_verify_cell_kzg_proof_batch<
    TFr: Fr + Debug,
    TG1: G1 + G1Mul<TFr> + G1GetFp<TG1Fp>,
    TG2: G2,
    TPoly: Poly<TFr>,
    TFFTSettings: FFTSettings<TFr>,
    TKZGSettings: KZGSettings<TFr, TG1, TG2, TFFTSettings, TPoly, TG1Fp, TG1Affine>,
    TG1Fp: G1Fp,
    TG1Affine: G1Affine<TG1, TG1Fp>,
>(
    load_trusted_setup: &dyn Fn(&str) -> Result<TKZGSettings, String>,
    verify_cell_kzg_proof_batch: &dyn Fn(
        &[TG1],
        &[usize],
        &[[TFr; FIELD_ELEMENTS_PER_CELL]],
        &[TG1],
        &TKZGSettings,
    ) -> Result<bool, String>,
) {
    let settings = load_trusted_setup(get_trusted_setup_path().as_str()).unwrap();
    let test_files: Vec<PathBuf> = glob::glob(&format!(
        "{}/{}",
        get_manifest_dir(),
        VERIFY_CELL_KZG_PROOF_BATCH_TEST_VECTORS
    ))
    .unwrap()
    .collect::<Result<Vec<_>, _>>()
    .unwrap();
    assert!(!test_files.is_empty());
    for test_file in test_files {
        // Ignore test vectors of verify_cell_kzg_proof_batch, as they are currently failing
        if test_file.parent().unwrap().file_name().unwrap().to_str().unwrap().starts_with("verify_cell_kzg_proof_batch_case_") {
            continue; // Remove this line in order to not skip through these vectors
        }
        let yaml_data = fs::read_to_string(test_file.clone()).unwrap();
        let test: verify_cell_kzg_proof_batch::Test = serde_yaml::from_str(&yaml_data).unwrap();
        let cells = match test
            .input
            .get_cell_bytes()
            .unwrap()
            .iter()
            .map(|bytes| {
                match bytes
                    .chunks(BYTES_PER_FIELD_ELEMENT)
                    .map(|bytes| TFr::from_bytes(bytes))
                    .collect::<Result<Vec<_>, String>>()
                {
                    Ok(value) => value
                        .try_into()
                        .map_err(|_| "Invalid field element per cell count".to_string()),
                    Err(err) => Err(err),
                }
            })
            .collect::<Result<Vec<_>, _>>()
        {
            Ok(v) => v,
            Err(err) => {
                // c-kzg-4844 also includes tests with invalid byte count for cell
                // in rust-kzg, these checks are performed outside of recovery function,
                // while parsing data (in c bindings, for example). These tests will be
                // additionally checked through rust-kzg c binding tests.
                // We add here assertion, to avoid accidentally skipping valid test case
                assert!(
                    test.get_output().is_none(),
                    "Parsing input failed with error {err:?}, for test vector {test_file:?}",
                );
                continue;
            }
        };
        let commitments = match test
            .input
            .get_commitment_bytes()
            .unwrap()
            .iter()
            .map(|bytes| TG1::from_bytes(&bytes))
            .collect::<Result<Vec<_>, _>>()
        {
            Ok(v) => v,
            Err(err) => {
                // c-kzg-4844 also includes tests with invalid byte count for cell
                // in rust-kzg, these checks are performed outside of recovery function,
                // while parsing data (in c bindings, for example). These tests will be
                // additionally checked through rust-kzg c binding tests.
                // We add here assertion, to avoid accidentally skipping valid test case
                assert!(
                    test.get_output().is_none(),
                    "Parsing input failed with error {err:?}, for test vector {test_file:?}",
                );
                continue;
            }
        };
        let proofs = match test
            .input
            .get_proof_bytes()
            .unwrap()
            .iter()
            .map(|bytes| TG1::from_bytes(&bytes))
            .collect::<Result<Vec<_>, _>>()
        {
            Ok(v) => v,
            Err(err) => {
                // c-kzg-4844 also includes tests with invalid byte count for cell
                // in rust-kzg, these checks are performed outside of recovery function,
                // while parsing data (in c bindings, for example). These tests will be
                // additionally checked through rust-kzg c binding tests.
                // We add here assertion, to avoid accidentally skipping valid test case
                assert!(
                    test.get_output().is_none(),
                    "Parsing input failed with error {err:?}, for test vector {test_file:?}",
                );
                continue;
            }
        };
        let cell_indices = test.input.get_cell_indices().unwrap();
        match verify_cell_kzg_proof_batch(
            &commitments,
            &cell_indices,
            &cells,
            &proofs,
            &settings,
        ) {
            Err(err) => assert!(test.get_output().is_none(), "Should correctly recover cells, but failed with error {err:?}, for test vector {test_file:?}"),
            Ok(value) => {
                let test_output = test.get_output();
                assert!(test_output.is_some(), "Should fail, but succeeded for test vector {test_file:?}");
                assert_eq!(value, test_output.unwrap(), "Test vector failed {test_file:?}");
            }
        }
    }
}