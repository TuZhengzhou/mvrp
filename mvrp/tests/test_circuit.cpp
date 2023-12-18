#include "structs.hpp"
#include "prover_circuit.hpp"
#include "verifier_circuit.hpp"
#include "circuit_generator.hpp"

using namespace cred;

void test_circuit() {

  libff::default_ec_pp::init_public_params();
    
    CredSRS srs;
    vector<vector<Fr>> WL;
    vector<vector<Fr>> WR;
    vector<vector<Fr>> WO;
    vector<vector<Fr>> WU;
    vector<Fr> c; 
    map<size_t, size_t> map;
    SET S;
    vector<Fr> aL;
    vector<Fr> aR;
    vector<Fr> aO;
    IPAProveSystem ipa_sys;

    size_t N = 1e4;
    size_t n, m;
    n = 32;
    m = 5;

    srs = CredSRS(N);
    
    auto gen = CircuitParaGenerator(n, m);
    
    WL = gen.W_L;
    WR = gen.W_R;
    WO = gen.W_O;
    WU = gen.W_V;
    c  = gen.c;
    S  = SET(gen.v);
    aL = gen.a_L;
    aR = gen.a_R;
    aO = gen.a_O;
    ipa_sys = IPAProveSystem();
    for(size_t i = 1; i <= m; i++) {
        map[i] = i;
    }

    ProofCircuit pi;
    Agenda agenda;

    ProverCircuit prover(srs, WL, WR, WO, WU, c, map, S, aL, aR, aO, ipa_sys);

    G1 V = prover._V;
    VerifierCircuit verifier(srs, V, WL, WR, WO, WU, c, map, ipa_sys);

    pi = prover.prove(agenda, false);
    bool result = verifier.verify(pi, agenda, false);
    // auto prover   = ProverCircuit(srs, set, ranges, ipa_sys);
    printf("result = %d\n", result);
}

int main() {
    test_circuit();
}