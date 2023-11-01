/** @file
*****************************************************************************

Declaration of Kate Polynomial Commitment in the generic group (GG) model.

This includes:
- class for commitment key: G1_vector
- class for commitment: G1 element
- class for KZGProof
- class for polynomial: Fr_vector
- PK generator algorithm
- commit algorithm
- create KZGProof/evaluation algorithm
- evaluation verifier algorithm

The implementation instantiates the protocol of \[KZG10].

Acronyms:

- vCNN+ = "Committed verifiable Convolutional Neural Network"

References:

\[KZG10]:
 "Polynomial Commitments",
 Aniket Kate, Gregory M. Zaverucha, Ian Goldberg
 ASIACRYPT 2010,
 <https://cacr.uwaterloo.ca/techreports/2010/cacr2010-10.pdf>

*****************************************************************************/
#include "structs.hpp"

#ifndef CRED_KZG_HPP_
#define CRED_KZG_HPP_

namespace cred
{

/***************************** Main algorithms *******************************/

/**
 * A setup algorithm for the KZG10.
 *
 * Given an authority t: degree, this algorithm produces commitment key, which is a t-SDH tuple.
 */


kzgCommitKey kzg_setup(const size_t t);

/**
 * Random Point Generator for the KZG10.
 *
 * Given three commitments, both Prover & Verifier generates random evaluation point: SHA256(Commit(A).x, Commit(B).x, Commit(C).x)
 * produces an hash = random evaluation point. of the polynomial. Non-interactively (Fiat-Shamir Heuristic)
 */


Fr kzg_hash(const G1 &commit_a, const G1 &commit_b, const G1 &commit_c);
Fr kzg_hash(const G1 &commit);

/**
 * A commit algorithm for the KZG10.
 *
 * Given a public key and polynomial, this algorithm
 * produces a commitment of the polynomial.
 */

// G1 kzg_commit(std::vector<G1> &ck, std::vector<Fr> &poly, size_t t);
G1 kzg_commit(const kzgCommitKey &ck, const std::vector<Fr> &poly, const size_t t);

/**
 * A KZGProof-generate algorithm for the KZG10.
 *
 * Given a public key, polynomial, and evaluation point, this algorithm produces a KZGProof of the evaluation of the polynomial.
 * (It proves that Polynomial is evaluated at particular evaluation point)
 */

KZGProof kzg_prove(const kzgCommitKey &ck, std::vector<Fr> &poly, const Fr &point, const size_t t);

pair<G1, KZGProof> kzg_prove(CredSRS& srs, const size_t exp_start, const size_t exp_span, std::vector<Fr>& poly, vector<G1>::iterator g1s_iter, vector<G2>::iterator g2s_iter);

/**
 * Polynomial Evaluation algorithm for the KZG10.
 *
 * Given a polynomial and point k, this algorithm evaluates the polynomial at point k.
 * "Polynomial is evaluated at particular evaluation point."
 */


Fr kzg_evaluate(const std::vector<Fr> &poly, const Fr &point, const size_t t);


 /**
 * A Evaluation Verifier algorithm for the KZG10.
 *
 * Given a public key, commitment,and KZGProof, this algorithm verifies the following statement.
 * "Polynomial is evaluated at particular evaluation point."
 */

bool kzg_vfyeval(const kzgCommitKey &ck, const G1 &commit, const KZGProof &KZGProof);

bool kzg_test();

} //cred

#endif // CRED_KZG_HPP_