#include "structs.hpp"

using namespace cred;

int main() {
  libff::start_profiling();
  libff::default_ec_pp::init_public_params();

  Fr a, b, c;
  G1 d, e, f;

  const int repeat = 1000;
  Agenda agenda;

  a = Fr::random_element();
  b = Fr::random_element();
  
  agenda.create_item("Fr_add");
  for(int i = 0; i < repeat; i++) {
    c = a + b;
  }
  agenda.mark_item_end("Fr_add");

  agenda.create_item("Fr_mul");
  for(int i = 0; i < repeat; i++) {
    c = a * b;
  }
  agenda.mark_item_end("Fr_mul");
  
  d = G1::random_element();
  e = G1::random_element();
  a = Fr::random_element();

  agenda.create_item("G1_mul");
  for(int i = 0; i < repeat; i++) {
    f = d + e;
  }
  agenda.mark_item_end("G1_mul");

  agenda.create_item("G1_exp");
  for(int i = 0; i < repeat; i++) {
    f = a * d;
  }
  agenda.mark_item_end("G1_exp");
  
  agenda.write_file("./microbenchmark.txt");

  return 0;
}
