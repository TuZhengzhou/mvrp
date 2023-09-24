/** @file
 *****************************************************************************

 Implementation of interfaces for Kate Commitment for R1CS example.

 *****************************************************************************
 * @author     This file is part of libsnark, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#ifndef RUN_KZG_TCC_
#define RUN_KZG_TCC_

#include <sstream>
#include <type_traits>

#include <libff/common/profiling.hpp>

#include <libsnark/zk_proof_systems/kzg/kzg.hpp>

namespace libsnark {

/**
 * The code below provides an example of all stages of running a Kate Commitment.
 *
 * This is the version we use for the DEMO CHECK
 */
template<typename ppT>
bool run_kzg(const r1cs_example<libff::Fr<ppT> > &example,
                        const bool test_serialization)
{

    libff::enter_block("Call to run_kzg");

    int t = 10000;

    /* Generate Polynomial to Commit: we need to put Convolution Poly. in this section */

    libff::Fr_vector<ppT> poly(t);

    for(int i = 0; i < t; i++) {
        libff::Fr<ppT> random = libff::Fr<ppT>::random_element();
        poly[i] = random;
        // poly[i].print();
    }

    /* Generate Random Point for Evaluation */

    libff::Fr<ppT> point = libff::Fr<ppT>::random_element();
    // libff::Fr<ppT> point = convert<ppT>(1);
    // point.print();

    /* Generate t-SDH tuple, and select secret randomness t */

    libff::print_header("Generate Key: t-SDH Tuple");
    commitkey<ppT> ck = kzg_setup<ppT>(t);
    printf("\n"); libff::print_indent(); libff::print_mem("after setup");

    /* Commit Polynomial into Product: G1-element */

    libff::print_header("Commit Polynomial");
    libff::G1<ppT> commit = kzg_commit<ppT>(ck, poly, t);
    printf("\n"); libff::print_indent(); libff::print_mem("after commit");

    /* Generate witness of the evaluation + Evaluate the Polynomial */

    libff::print_header("Create Witness");
    witness<ppT> wit = kzg_witness<ppT>(ck, poly, point, t);
    printf("\n"); libff::print_indent(); libff::print_mem("after create-witness");

    /* Verify evaluation */
    libff::print_header("Verify Evaluation of Polynomial");
    bool verifyresult = kzg_vfyeval<ppT>(ck, commit, wit);

    if (verifyresult == true) {
        libff::print_header("VERIFICATION ACCEPT!!");
    } else {
        libff::print_header("VERIFICATION REJECT");
    }
    
    printf("\n"); libff::print_indent(); libff::print_mem("after vfyeval");

    /* Verify evaluation TEST-DEMO*/
    // libff::print_header("Verify TEST-Logical Error");
    // bool testresult = kzg_testvfy<ppT>(ck, commit, wit, poly2, t);
    // printf("\n"); libff::print_indent(); libff::print_mem("after testvfy");

    //clear all
    poly.clear();
    point.clear();

    return verifyresult;

    libff::leave_block("Call to run_kzg");
    
}

} // libsnark

#endif // RUN_KZG_TCC_