#include "structs.hpp"
// #include "basic_types.hpp"
#include "circuit_prover.hpp"
#include "circuit_verifier.hpp"
#include "circuit_generator.hpp"
#include "linear_combination.hpp"

using namespace cred;

bool test_circuit() {

  libff::default_ec_pp::init_public_params();

  CCredSRS srs;
  vector<LinearCombination> constraints;
  map<size_t, size_t> map;
  CSet S;
  vector<Fr> aL;
  vector<Fr> aR;
  vector<Fr> aO;
  IPAProveSystem ipa_sys;

  size_t N = (int)pow(2,6) + 1;
  size_t num_multi_gates, num_commit_input;
  num_multi_gates = (int)pow(2,6);
  num_commit_input = 5;

  srs = CCredSRS(N);
  
  auto gen = CircuitParaGenerator(num_multi_gates, num_commit_input);
  
  constraints = gen.constraints;
  S  = CSet(gen.v);
  aL = gen.a_L;
  aR = gen.a_R;
  aO = gen.a_O;
  ipa_sys = IPAProveSystem();
  for(size_t i = 1; i <= num_commit_input; i++) {
      map[i] = i;
  }

  CCircuitProof pi;
  Agenda agenda;

  ProverCircuit prover(srs, constraints, num_multi_gates, num_commit_input, map, S, aL, aR, aO, ipa_sys);

  G1 V = prover._V;
  VerifierCircuit verifier(srs, constraints, V, num_multi_gates, num_commit_input, map, ipa_sys);

  libff::enter_block("ProverCircuit");
  pi = prover.prove(agenda, false);
  libff::leave_block("ProverCircuit");

  libff::enter_block("VerifierCircuit");
  // bool result = verifier.verify(pi, agenda, false);
  bool result = verifier.verify(pi, agenda, true);
  libff::leave_block("VerifierCircuit");

  // auto prover   = ProverCircuit(srs, set, ranges, ipa_sys);
  printf("result = %d\n", result);

  libff::print_cumulative_times();
  libff::print_mem("after all");
  // libff::print_cumulative_time_entry("ProverCircuit");
  // libff::print_cumulative_time_entry("VerifierCircuit");

  return result;
}

int main() {
  bool result = test_circuit();
  return (int)(!result);
}