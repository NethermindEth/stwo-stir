//! Implements a FRI polynomial commitment scheme.
//! This is a protocol where the prover can commit on a set of polynomials and then prove their
//! opening on a set of points.
//! Note: This implementation is not really a polynomial commitment scheme, because we are not in
//! the unique decoding regime. This is enough for a STARK proof though, where we only want to imply
//! the existence of such polynomials, and are ok with having a small decoding list.
//! Note: Opened points cannot come from the commitment domain.

mod prover;
pub mod quotients;
mod utils;
mod verifier;

use crate::core::backend::Backend;
use crate::core::channel::Blake2sChannel;
use crate::core::circle::CirclePoint;
use crate::core::ColumnVec;
use crate::core::fields::qm31::SecureField;
use crate::core::fri::{FriConfig, FriOps, FriProof};
use crate::core::poly::circle::PolyOps;
use crate::core::prover::{VerificationError, LOG_BLOWUP_FACTOR, LOG_LAST_LAYER_DEGREE_BOUND, N_QUERIES};
use crate::core::vcs::blake2_merkle::Blake2sMerkleHasher;
use crate::core::vcs::ops::MerkleOps;
pub use self::prover::{CommitmentSchemeProof, CommitmentSchemeProver, CommitmentTreeProver};
pub use self::utils::TreeVec;
pub use self::verifier::CommitmentSchemeVerifier;

#[derive(Copy, Debug, Clone)]
pub struct TreeColumnSpan {
    pub tree_index: usize,
    pub col_start: usize,
    pub col_end: usize,
}

pub trait PolynomialCommitmentSchemeBase {
    type Config;
    type Proof;
    type Channel;

    fn config() -> Self::Config;
}

pub trait PolynomialCommitmentScheme<B: PolyOps>: PolynomialCommitmentSchemeBase {
    type Prover<'a>: PolynomialProver<Self::Channel, Self::Proof> where B: 'a;
    type Verifier: PolynomialVerifier<Self::Channel, Self::Proof>;
}

pub trait PolynomialProver<Channel, Proof> {
    fn prove_values(
        &self,
        sampled_points: TreeVec<ColumnVec<Vec<CirclePoint<SecureField>>>>,
        channel: &mut Channel,
    ) -> CommitmentSchemeProof<Proof>;
}

pub trait PolynomialVerifier<Channel, Proof>: Sized {
    fn verify_values(
        &self,
        sampled_points: TreeVec<ColumnVec<Vec<CirclePoint<SecureField>>>>,
        proof: CommitmentSchemeProof<Proof>,
        channel: &mut Channel,
    ) -> Result<(), VerificationError>;
}

// TODO(alex): Move to a submodule
pub struct FriCommitmentScheme;

impl PolynomialCommitmentSchemeBase for FriCommitmentScheme {
    type Config = FriConfig;
    type Proof = FriProof<Blake2sMerkleHasher>;
    type Channel = Blake2sChannel;

    fn config() -> FriConfig {
        FriConfig::new(LOG_LAST_LAYER_DEGREE_BOUND, LOG_BLOWUP_FACTOR, N_QUERIES)
    }
}

impl<B> PolynomialCommitmentScheme<B> for FriCommitmentScheme
where
    B: Backend + FriOps + MerkleOps<Blake2sMerkleHasher>,
{
    type Prover<'a> = CommitmentSchemeProver<'a, B, Self::Proof> where B: 'a;
    type Verifier = CommitmentSchemeVerifier<Self::Proof>;
}
