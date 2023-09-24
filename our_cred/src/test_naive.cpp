#include <iostream>
#include <vector>
#include <map>
#include <assert.h>
#include <iterator>
#include "libff/common/profiling.hpp"
#include "libff/algebra/curves/mnt/mnt6/mnt6_pp.hpp"
#include "libff/algebra/curves/bn128/bn128_pp.hpp"
#include "libff/algebra/curves/public_params.hpp"
#include "structs.hpp"
#include "prover.hpp"
#include "verifier.hpp"
#include "ipa.hpp"

/*
  
*/

using namespace std;
using namespace cred;
using cred::G1;
using cred::G2;
using cred::GT;
using cred::Fr;

int main() {
    // libff::start_profiling();
    libff::default_ec_pp::init_public_params();

    size_t N = (size_t)1 << 20;
    size_t max_bit_len = (size_t)1 << 8;
    CredSRS srs = CredSRS(N);

    size_t set_size = max((size_t)5, N);
    std::vector<Fr> set_values_vec;
    for(size_t i = 0; i < set_size; i++) {
      set_values_vec.emplace_back(Fr::one() + Fr::one());
    }
    SET set = SET(set_values_vec);

    std::vector<size_t> I = {1,2,3,4};
    std::vector<Range> range_vec;

    for(size_t i = 0; i < I.size(); i++) {
        // range_vec.emplace_back(Range(Fr::zero(), -Fr::one(), Fr::num_bits));
        range_vec.emplace_back(Range(Fr::zero(), -Fr::one(), size_t(4)));
    }
    Ranges ranges = Ranges(I, range_vec);

    G1 u = Fr::random_element() * G1::one();
    IPAProveSystem ipa_sys = IPAProveSystem(u);

    Prover prover = Prover(srs, set, ranges, ipa_sys);
    CredProof pi = prover.prove();
    G1 C = prover._C;

    Verifier verifier = Verifier(srs, C, set_size, ranges, ipa_sys);
    bool result = verifier.verify(pi, prover);
    cout << "verifier result: " << result << endl;
    cout << "size of pi: " << sizeof(pi) << endl;
    cout << "size of pi_tilde: " << sizeof(pi._pi_tilde) << endl;


    // IPAOneRecursionTest();
    // IPATest();
    // IPATest(64);

    cout << "output OK:>" << endl;
}
