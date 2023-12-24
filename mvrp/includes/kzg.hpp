/** @file
*****************************************************************************

Declaration of Kate Polynomial Commitment in the generic group (GG) model.

This includes:
- class for commitment key: G1_vector
- class for commitment: G1 element
- class for CKzgProof
- class for polynomial: Fr_vector
- PK generator algorithm
- commit algorithm
- create CKzgProof/evaluation algorithm
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

#ifndef CRED_KZG_HPP_
#define CRED_KZG_HPP_

#include <vector>
#include <utility>
#include <cstddef>
#include <xassert/XAssert.h>
#include "basic_types.hpp"

using std::vector;
using std::pair;

namespace cred
{

/******************************** Commitment key ********************************/
// class CKzgKey;
// std::ostream& operator<<(std::ostream &out, const CKzgKey &ck);
// std::istream& operator>>(std::istream &in, CKzgKey &ck);

class CKzgKey
{
    public:
    std::vector<G1> g1;
    std::vector<G2> g2;

    CKzgKey() = default;
    CKzgKey& operator=(const CKzgKey &other) = default;
    CKzgKey(const CKzgKey &other) = default;
    CKzgKey(CKzgKey &&other) = default;
    CKzgKey(
        std::vector<G1> &&g1,
        std::vector<G2> &&g2) :
    g1(std::move(g1)),
    g2(std::move(g2))
    {};
    CKzgKey(
        std::vector<G1> &g1,
        std::vector<G2> &g2) :
    g1(g1),
    g2(g2)
    {};

    size_t G1_size() const
    {
        return g1.size();
    }

    size_t G2_size() const
    {
        return g2.size();
    }

    size_t GT_size() const
    {
        return 1;
    }

    size_t size_in_bits() const
    {
      size_t ret;
      ret += g1.size() > 0 ? g1.size() * sizeof(g1[0]) : 0;
      ret += g2.size() > 0 ? g2.size() * sizeof(g2[0]) : 0;
      return ret;
    }

    void print_size() const
    {
        libff::print_indent(); printf("* G1 elements in CommitKey: %zu\n", this->G1_size());
        libff::print_indent(); printf("* G2 elements in CommitKey: %zu\n", this->G2_size());
        libff::print_indent(); printf("* Commit Key size in bits: %zu\n", this->size_in_bits());
    }

    bool operator==(const CKzgKey &other) const;
    friend std::ostream& operator<< (std::ostream &out, const CKzgKey &ck);
    friend std::istream& operator>> (std::istream &in, CKzgKey &ck);
};

/******************************** Witness ********************************/

class CKzgProof
{
    public:
    // Fr point; // can be computed
    G1 eval;    // [f(point)]_1
    std::vector<Fr> psi;
    G1 w1;      // [h(point) except the N+1 term]_1
    GT wt;      // [the N+1 term of h(point)]_2, set to 0 when h(x) has no N+1 term

    CKzgProof() = default;
    CKzgProof& operator=(const CKzgProof &other) = default;
    CKzgProof(const CKzgProof &other) = default;
    CKzgProof(CKzgProof &&other) = default;
    CKzgProof(
        G1 &&eval,
        G1 &&w1):

    // point(std::move(point)),
    eval(std::move(eval)),
    w1(std::move(w1))
    {};

    size_t G1_size() const
    {
        return sizeof(w1) + sizeof(eval);
    }

    size_t size_in_bits() const
    {
       return (1 + w1.size_in_bits() + eval.size_in_bits());
    }

    void print_size() const
    {
        libff::print_indent(); printf("* G1 elements in Witness: %zu\n", this->G1_size());
        libff::print_indent(); printf("* Witness size in bits: %zu\n", this->size_in_bits());
    }

    // bool operator==(const CKzgProof &other) const;
    // friend std::ostream& operator<< (std::ostream &out, const CKzgProof &wit);
    // friend std::istream& operator>> (std::istream &in, CKzgProof &wit);
};

// std::ostream& operator<< (std::ostream &out, const CKzgProof &wit);
// std::istream& operator>> (std::istream &in, CKzgProof &wit);

/***************************** Main algorithms *******************************/

/**
 * A setup algorithm for the KZG10.
 *
 * Given an authority t: degree, this algorithm produces commitment key, which is a t-SDH tuple.
 */


CKzgKey kzg_setup(const size_t t);

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
G1 kzg_commit(const CKzgKey &ck, const std::vector<Fr> &poly, const size_t t);

/**
 * A CKzgProof-generate algorithm for the KZG10.
 *
 * Given a public key, polynomial, and evaluation point, this algorithm produces a CKzgProof of the evaluation of the polynomial.
 * (It proves that Polynomial is evaluated at particular evaluation point)
 */

CKzgProof kzg_prove(const CKzgKey &ck, std::vector<Fr> &poly, const Fr &point, size_t t);

CKzgProof kzg_prove(const CCredSRS& srs, std::vector<Fr>& poly, size_t t, G1& commit);

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
 * Given a public key, commitment,and CKzgProof, this algorithm verifies the following statement.
 * "Polynomial is evaluated at particular evaluation point."
 */

bool kzg_vfyeval(const CKzgKey &ck, const G1 &commit, const CKzgProof &CKzgProof);

bool kzg_test();

} // namespace cred

#endif // CRED_KZG_HPP_