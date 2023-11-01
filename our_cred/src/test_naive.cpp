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
#include "kzg.hpp"
#include "circuit_generator.hpp"

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

    auto gen = CircuitParaGenerator(50, 40);
    bool check = gen.check();

    cout << "output OK:>" << endl;
}
