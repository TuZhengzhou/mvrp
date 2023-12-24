#include "libff/common/default_types/ec_pp.hpp"
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

    return (int)(!check);
}
