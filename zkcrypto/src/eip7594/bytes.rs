use crate::kzg_types::{Fr, G1, C_KZG_RET}; 
use kzg::eip_4844::{C_KZG_RET_BADARGS, C_KZG_RET_OK, Bytes32};
use zkcrypto::lib::Scalar;
use blst::{blst_fr_from_scalar, blst_scalar_fr_check, blst_scalar_from_bendian};

pub fn bytes_to_bls_field(out: &mut Fr, b: &Bytes32) -> C_KZG_RET {
    let mut tmp = Scalar::default();
    blst_scalar_from_bendian(&mut tmp, &b.bytes);
    if !blst_scalar_fr_check(&tmp) {
        return C_KZG_RET_BADARGS;
    }
    blst_fr_from_scalar(out, &tmp);
    return C_KZG_RET_OK
}
