#ifndef CRED_IPA_HPP
#define CRED_IPA_HPP
#include <vector>
#include <map>
#include <tuple>
#include <xassert/XAssert.h>
#include "libff/common/profiling.hpp"
#include "libff/algebra/curves/public_params.hpp"
#include "libff/common/default_types/ec_pp.hpp"
#include "structs.hpp"


using namespace std;

namespace cred {

class IPAProveSystem {
public:

  G1 _u;
  G1 _P;
  Fr _c;

  IPAProveSystem() {};
  IPAProveSystem(G1& u);
  IPAProveSystem(G1& P, G1& u, Fr& c);
  G1 IPACommit();
  // IPAProveSystem() {};


  IPAProofOneRecursion IpaProveOneRecursion(
    const std::vector<G1>& g_vec, const std::vector<G1>& h_vec, \
    const std::vector<Fr>& a_vec, const std::vector<Fr>& b_vec, \
    const Fr& c
  );

  bool IpaVerifyOneRecursion(
    const IPAProofOneRecursion& pi, const std::vector<G1>& g_vec, const std::vector<G1>& h_vec
  );


  IPAProof IpaProve(const std::vector<G1>& g_vec, const std::vector<G1>& h_vec, \
                    const std::vector<Fr>& a_vec, const std::vector<Fr>& b_vec, \
                    const Fr& c);

  bool IpaVerify(const IPAProof& pi, const std::vector<G1>& g_vec, const std::vector<G1>& h_vec);

  bool IpaMultiExpVerify(const IPAProof& pi, const std::vector<G1>& g_vec, const std::vector<G1>& h_vec);

private:
  Fr generate_random(const Fr& pre_random, const G1& L, const G1& R);

  G1 ipa_hash(
    const std::vector<G1>& g_vec, const std::vector<G1>& h_vec, \
    const std::vector<Fr>& al_vec, const std::vector<Fr>& ar_vec, \
    const std::vector<Fr>& bl_vec, const std::vector<Fr>& br_vec, \
    const Fr& c
  );

  G1 ipa_hash(
    const std::vector<G1>& g_vec, const std::vector<G1>& h_vec, \
    const Fr& al, const Fr& ar, const Fr& bl, const Fr& br, const Fr& c
  );
  
  G1 ipa_hash(
    std::vector<G1>::const_iterator g_iter , std::vector<G1>::const_iterator h_iter, \
    std::vector<Fr>::const_iterator al_iter, std::vector<Fr>::const_iterator ar_iter, \
    std::vector<Fr>::const_iterator bl_iter, std::vector<Fr>::const_iterator br_iter, \
    const Fr& c, const size_t n_apos
  );
};

bool IPAOneRecursionTest();
bool IPATest(size_t n = 2);
};
#endif