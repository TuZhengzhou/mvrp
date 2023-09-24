/** @file
 *****************************************************************************

 Declaration of functionality that runs the Kate Polynomial Commitment for R1CS example.

 *****************************************************************************
 * @author     This file is part of libsnark, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#ifndef RUN_KZG_HPP_
#define RUN_KZG_HPP_

#include <libff/algebra/curves/public_params.hpp>

#include <libsnark/relations/constraint_satisfaction_problems/r1cs/examples/r1cs_examples.hpp>

namespace libsnark {

/**
 * Runs the Polynomial Commitment Scheme for a given
 * R1CS example (specified by a constraint system, input, and witness).
 *
 * Optionally, also test the serialization routines for keys and commitments.
 * (This takes additional time.)
 * I'm not actually sure should we do that. lol
 */
template<typename ppT>
bool run_kzg(const r1cs_example<libff::Fr<ppT> > &example,
                        const bool test_serialization);

} // libsnark

#include <libsnark/zk_proof_systems/kzg/examples/run_kzg.tcc>

#endif // RUN_KZG_HPP_
