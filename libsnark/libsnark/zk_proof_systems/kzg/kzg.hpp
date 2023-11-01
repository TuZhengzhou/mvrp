/** @file
*****************************************************************************

Declaration of Kate Polynomial Commitment in the generic group (GG) model.

This includes:
- class for commitment key: G1_vector
- class for commitment: G1 element
- class for witness
- class for polynomial: Fr_vector
- PK generator algorithm
- commit algorithm
- create witness/evaluation algorithm
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

*****************************************************************************
* @author     This file is part of libsnark, developed by SCIPR Lab
*             and contributors (see AUTHORS).
* @copyright  MIT license (see LICENSE file)
*****************************************************************************/

#ifndef KZG_HPP_
#define KZG_HPP_

namespace libsnark
{

/******************************** Commitment key ********************************/
template<typename ppT>
class commitkey;

template<typename ppT>
std::ostream& operator<<(std::ostream &out, const commitkey<ppT> &ck);

template<typename ppT>
std::istream& operator>>(std::istream &in, commitkey<ppT> &ck);

template <typename ppT>
class commitkey
{
    public:
    libff::G1_vector<ppT> g1;
    libff::G2_vector<ppT> g2;

    commitkey() = default;
    commitkey<ppT>& operator=(const commitkey<ppT> &other) = default;
    commitkey(const commitkey<ppT> &other) = default;
    commitkey(commitkey<ppT> &&other) = default;
    commitkey(
        libff::G1_vector<ppT> &&g1,
        libff::G2_vector<ppT> &&g2) :
    g1(std::move(g1)),
    g2(std::move(g2))
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
       return (g1.size_in_bits() + g2.size_in_bits());
    }

    void print_size() const
    {
        libff::print_indent(); printf("* G1 elements in CommitKey: %zu\n", this->G1_size());
        libff::print_indent(); printf("* G2 elements in CommitKey: %zu\n", this->G2_size());
        libff::print_indent(); printf("* Commit Key size in bits: %zu\n", this->size_in_bits());
    }

    bool operator==(const commitkey<ppT> &other) const;
    friend std::ostream& operator<< <ppT>(std::ostream &out, const commitkey<ppT> &ck);
    friend std::istream& operator>> <ppT>(std::istream &in, commitkey<ppT> &ck);
};

/******************************** Witness ********************************/
template<typename ppT>
class witness;

template<typename ppT>
std::ostream& operator<<(std::ostream &out, const witness<ppT> &wit);

template<typename ppT>
std::istream& operator>>(std::istream &in, witness<ppT> &wit);

template <typename ppT>
class witness
{
    public:
    libff::Fr<ppT> point;
    libff::G1<ppT> eval;
    // libff::Fr_vector<ppT> psi; //need to be deleted
    libff::G1<ppT> w1;
    // libff::G2<ppT> w2; //need to be deleted.

    witness() = default;
    witness<ppT>& operator=(const witness<ppT> &other) = default;
    witness(const witness<ppT> &other) = default;
    witness(witness<ppT> &&other) = default;
    witness(
        libff::Fr<ppT> &&point,
        libff::G1<ppT> &&eval,
        libff::G1<ppT> &&w1):

    point(std::move(point)),
    eval(std::move(eval)),
    w1(std::move(w1))
    {};

    size_t G1_size() const
    {
        return w1.size() + eval.size();
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

    bool operator==(const witness<ppT> &other) const;
    friend std::ostream& operator<< <ppT>(std::ostream &out, const witness<ppT> &wit);
    friend std::istream& operator>> <ppT>(std::istream &in, witness<ppT> &wit);
};

/***************************** Main algorithms *******************************/

/**
 * A setup algorithm for the KZG10.
 *
 * Given an authority t: degree, this algorithm produces commitment key, which is a t-SDH tuple.
 */

template<typename ppT>
commitkey<ppT> kzg_setup(int t);

/**
 * Random Point Generator for the KZG10.
 *
 * Given three commitments, both Prover & Verifier generates random evaluation point: SHA256(Commit(A).x, Commit(B).x, Commit(C).x)
 * produces an hash = random evaluation point. of the polynomial. Non-interactively (Fiat-Shamir Heuristic)
 */

template<typename ppT>
libff::Fr<ppT> kzg_hash(libff::G1<ppT> &commit_a, libff::G1<ppT> &commit_b, libff::G1<ppT> &commit_c);

/**
 * A commit algorithm for the KZG10.
 *
 * Given a public key and polynomial, this algorithm
 * produces a commitment of the polynomial.
 */
template<typename ppT>
// libff::G1<ppT> kzg_commit(libff::G1_vector<ppT> &ck, libff::Fr_vector<ppT> &poly, int t);
libff::G1<ppT> kzg_commit(commitkey<ppT> &ck, libff::Fr_vector<ppT> &poly, int t);

/**
 * A witness-generate algorithm for the KZG10.
 *
 * Given a public key, polynomial, and evaluation point, this algorithm produces a witness of the evaluation of the polynomial.
 * (It proves that Polynomial is evaluated at particular evaluation point)
 */
template<typename ppT>
witness<ppT> kzg_witness(commitkey<ppT> &ck, libff::Fr_vector<ppT> &poly, libff::Fr<ppT> &point, int t);

/**
 * Polynomial Evaluation algorithm for the KZG10.
 *
 * Given a polynomial and point k, this algorithm evaluates the polynomial at point k.
 * "Polynomial is evaluated at particular evaluation point."
 */

template<typename ppT>
libff::Fr<ppT> kzg_evaluate(libff::Fr_vector<ppT> &poly, libff::Fr<ppT> &point, int t);


 /**
 * A Evaluation Verifier algorithm for the KZG10.
 *
 * Given a public key, commitment,and witness, this algorithm verifies the following statement.
 * "Polynomial is evaluated at particular evaluation point."
 */
template<typename ppT>
bool kzg_vfyeval(commitkey<ppT> &ck, libff::G1<ppT> &commit, witness<ppT> &witness);


} //libsnark

#include <libsnark/zk_proof_systems/kzg/kzg.tcc>

#endif // KZG_HPP_