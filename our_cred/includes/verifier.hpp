#ifndef CRED_VERIFIER_HPP
#define CRED_VERIFIER_HPP
#include <vector>
#include <map>
#include <tuple>
#include "libff/algebra/curves/public_params.hpp"
#include "libff/algebra/fields/field_utils.hpp"
#include "structs.hpp"

namespace cred {

class Verifier{
public:
  CredSRS _srs; // srs
  size_t _set_size;   // set S
  G1 _C;    // commitment C
  Ranges _ranges;
  IPAProveSystem _ipa_prove_sys;

  size_t _D;    // D = sum_(j in I) n_(j)
  size_t _N;    // _N = srs.N
  map<size_t, size_t> _part_sums; // sum_(i=1)^(j-1) n_i

  Fr _y, _z, _x, _s;  // randoms used in interaction
  
  Fr _delta;
  GT _left_1, _left_2, _left_3, _left;
  GT _right_1, _right_2, _right_3, _right;

  CredProof _pi_cred;

  Verifier() {};
  Verifier(const CredSRS srs, const G1 C, const size_t set_size, const Ranges& ranges, const IPAProveSystem& ipa_sys);

  bool verify(const CredProof& pi, const Prover& prover);

private:

};

  
};

#endif