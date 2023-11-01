/** @file
 *****************************************************************************
 Test program that exercises the KZG10 (first generator, then
 prover, then verifier) on a synthetic R1CS instance.

 *****************************************************************************
 * @author     This file is part of libsnark, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/
#include <cassert>
#include <cstdio>

#include <libff/common/profiling.hpp>
#include <libff/common/utils.hpp>

#include <libsnark/common/default_types/kzg_pp.hpp>
#include <libsnark/relations/constraint_satisfaction_problems/r1cs/examples/r1cs_examples.hpp>
#include <libsnark/zk_proof_systems/kzg/examples/run_kzg.hpp>

using namespace libsnark;

template<typename ppT>
void test_kzg(size_t num_constraints,
              size_t input_size)
{
    libff::print_header("(enter) Test Kate Commitment");

    const bool test_serialization = true;
    // r1cs_example<libff::Fr<ppT> > example = generate_r1cs_example_with_binary_input<libff::Fr<ppT> >(num_constraints, input_size);
    r1cs_example<libff::Fr<ppT> > example = generate_r1cs_example_with_binary_input<libff::Fr<ppT> >(1, 1);
    const bool bit = run_kzg<ppT>(example, test_serialization);
    assert(bit);

    libff::print_header("(leave) Test Kate Commitment");
}

int main()
{
    default_kzg_pp::init_public_params();
    libff::start_profiling();

    test_kzg<default_kzg_pp>(1000, 100);
}
