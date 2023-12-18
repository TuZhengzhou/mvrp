#include "structs.hpp"

using namespace cred;

int main() {
  libff::start_profiling();
  libff::default_ec_pp::init_public_params();
  libff::inhibit_profiling_info = true;
  // libff::inhibit_profiling_counters = true;

  Fr a, b, c;
  G1 d, e, f;

  const int repeat = 1e4;

  a = Fr::random_element();
  b = Fr::random_element();
  
  libff::enter_block("Fr_add");
  for(int i = 0; i < repeat; i++) {
    c = a + b;
  }
  libff::leave_block("Fr_add");

  libff::enter_block("Fr_mul");
  for(int i = 0; i < repeat; i++) {
    c = a * b;
  }
  libff::leave_block("Fr_mul");
  
  d = G1::random_element();
  e = G1::random_element();
  a = Fr::random_element();

  libff::enter_block("G1_mul");
  for(int i = 0; i < repeat; i++) {
    f = d + e;
  }
  libff::leave_block("G1_mul");

  libff::enter_block("G1_exp");
  for(int i = 0; i < repeat; i++) {
    f = a * d;
  }
  libff::leave_block("G1_exp");

  libff::print_compilation_info();
  libff::print_cumulative_times();
  libff::print_cumulative_op_counts();
  libff::print_cumulative_op_counts(true);


  return 0;
}
