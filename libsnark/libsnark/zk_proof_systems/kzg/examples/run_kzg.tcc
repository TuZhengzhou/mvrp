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
#include <libff/algebra/fields/field_utils.hpp>
#include <libff/common/utils.hpp>
#include <libsnark/gadgetlib1/gadgets/hashes/sha256/sha256_gadget.hpp>

#include <libsnark/zk_proof_systems/kzg/kzg.hpp>

namespace libsnark {

/**
 * The code below provides an example of all stages of running a Kate Commitment.
 */
template<typename ppT>
bool run_kzg(const r1cs_example<libff::Fr<ppT> > &example,
                        const bool test_serialization)
{

    libff::enter_block("Call to run_kzg");

    /* Degree of Polynomial A, B, C */

    int t_a = 100;
    int t_b = 100;
    int t_c = t_a + t_b - 1;

    /* Generate Polynomial to Commit: we need to put Convolution Poly. in this section */

    libff::Fr_vector<ppT> poly_a(t_a);
    libff::Fr_vector<ppT> poly_b(t_b);
    libff::Fr_vector<ppT> poly_c(t_c);

    for(int i = 0; i < t_a; i++) {
        libff::Fr<ppT> random = libff::Fr<ppT>::random_element();
        poly_a[i] = random;
    }

    for(int i = 0; i < t_b; i++) {
        libff::Fr<ppT> random = libff::Fr<ppT>::random_element();
        poly_b[i] = random;
    }

    // C(x) = A(x) * B(x)
    for(int i = 0; i < t_a; i++) {
        for(int j = 0; j < t_b; j++) {
            poly_c[i + j] += poly_a[i] * poly_b[j];
        }
    }

    /* Generate t-SDH tuple, and select secret randomness t */

    libff::print_header("Generate Key: t-SDH Tuple");
    commitkey<ppT> ck = kzg_setup<ppT>(t_c);
    printf("\n"); libff::print_indent(); libff::print_mem("after setup");

    /* Commit Polynomial into Product: G1-element */

    libff::print_header("Commit Polynomial: A(x)");
    libff::G1<ppT> commit_a = kzg_commit<ppT>(ck, poly_a, t_a);
    printf("\n"); libff::print_indent(); libff::print_mem("after commit");

    libff::print_header("Commit Polynomial: B(x)");
    libff::G1<ppT> commit_b = kzg_commit<ppT>(ck, poly_b, t_b);
    printf("\n"); libff::print_indent(); libff::print_mem("after commit");

    libff::print_header("Commit Polynomial: C(x)");
    libff::G1<ppT> commit_c = kzg_commit<ppT>(ck, poly_c, t_c);
    printf("\n"); libff::print_indent(); libff::print_mem("after commit");

    /* Generate Random Point for Evaluation: SHA256 */

    libff::print_header("Generate Random point: Hash(Commit(A(x)), Commit(B(x)), Commit(C(x)))");
    libff::Fr<ppT> point = kzg_hash<ppT>(commit_a, commit_b, commit_c);
    printf("\n"); libff::print_indent(); libff::print_mem("after hash");

    /* Generate Evaluation Value for Convolution Proof: Groth16 */

    libff::Fr<ppT> eval_a = kzg_evaluate<ppT>(poly_a, point, t_a);
    libff::Fr<ppT> eval_b = kzg_evaluate<ppT>(poly_b, point, t_b); 
    libff::Fr<ppT> eval_c = kzg_evaluate<ppT>(poly_c, point, t_c);

    /* Generate witness of the evaluation + Evaluate the Polynomial */

    libff::print_header("Create Witness: A(x)");
    witness<ppT> wit_a = kzg_witness<ppT>(ck, poly_a, point, t_a);
    printf("\n"); libff::print_indent(); libff::print_mem("after create-witness");

    libff::print_header("Create Witness: B(x)");
    witness<ppT> wit_b = kzg_witness<ppT>(ck, poly_b, point, t_b);
    printf("\n"); libff::print_indent(); libff::print_mem("after create-witness");

    libff::print_header("Create Witness: C(x)");
    witness<ppT> wit_c = kzg_witness<ppT>(ck, poly_c, point, t_c);
    printf("\n"); libff::print_indent(); libff::print_mem("after create-witness");


    /* Verify evaluation */
    bool verifyresult;

    libff::print_header("Verify Evaluation of Polynomial: A(x)");
    bool verifyresult_a = kzg_vfyeval<ppT>(ck, commit_a, wit_a);

    libff::print_header("Verify Evaluation of Polynomial: B(x)");
    bool verifyresult_b = kzg_vfyeval<ppT>(ck, commit_b, wit_b);

    libff::print_header("Verify Evaluation of Polynomial: C(x)");
    bool verifyresult_c = kzg_vfyeval<ppT>(ck, commit_c, wit_c);

    if (verifyresult_a == true && verifyresult_b == true && verifyresult_c == true) {
        verifyresult = true;
        libff::print_header("VERIFICATION ACCEPT!!");
    } else {
        verifyresult = false;
        libff::print_header("VERIFICATION REJECT");
    }
    
    printf("\n"); libff::print_indent(); libff::print_mem("after vfyeval");

    /* Output Values (a(k), b(k), c(k)) to link with Convolution Proof: Groth16 */
    libff::print_header("vCNN+: Convolution Proof");

    libff::Fr_vector<ppT> jsnark(3);
    jsnark[0] = eval_a; jsnark[1] = eval_b; jsnark[2] = eval_c;
    printf("This values are for jsnark\n");
    jsnark[0].print(); jsnark[1].print(); jsnark[2].print();

    //clear all
    poly_a.clear();
    poly_b.clear();
    poly_c.clear();
    jsnark.clear();

    return verifyresult;

    libff::leave_block("Call to run_kzg");
    
}

} // libsnark

#endif // RUN_KZG_TCC_