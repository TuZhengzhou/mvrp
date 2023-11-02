#include "setup.hpp"
#include "structs.hpp"

namespace cred {
CredSRS trustSetup(uint64_t N) {
  assert(N >= 1);
  Fr beta = Fr::random_element();
  G1 g1_one = G1::random_element();
  G2 g2_one = G2::random_element();
  GT gt_one = ReducedPairing(g1_one, g2_one);

  CredSRS ret = CredSRS();
  ret.N = N;

  Fr exp = beta;
  for(int i = 1; i <= N; i++) {
    ret.g1s.emplace_back(exp * g1_one);
    ret.g2s.emplace_back(exp * g2_one);
    exp *= beta;
  }
  ret.gt = gt_one ^ exp;  // now exp is beta^(N+1)

  exp *= beta;    // get beta^(N+2)
  for(int i = N+2; i <= 2*N; i++) {
    ret.g1s.emplace_back(exp * g1_one);
    exp *= beta;
  }

  assert(ret.g1s.size() == 2*N-1);
  assert(ret.g2s.size() == N);
  assert(ret.gt = ret.g1s[N-1] * ret.g2s[0]);
  
  return ret;
}
};