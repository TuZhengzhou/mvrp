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
  CCredSRS _srs; // srs
  size_t _set_size;   // set S
  G1 _C;    // commitment C
  CRanges _ranges;
  IPAProveSystem _ipa_prove_sys;
  
  size_t _D;    // D = sum_(j in I) n_(j)
  size_t _N;    // _N = srs.N
  map<size_t, size_t> _part_sums; // sum_(i=1)^(j-1) n_i

  // CRangeProof _pi_cred;

  Verifier() {};
  Verifier(const CCredSRS& srs, const G1& C, const size_t& set_size, const CRanges& ranges, const IPAProveSystem& ipa_sys);

  bool verify(const CRangeProof& pi, Agenda& agenda, const bool improved = true);
  bool verify_base(const CRangeProof& pi, Agenda& agenda);
  bool verify_improved(const CRangeProof& pi, Agenda& agenda);

private:

};

  
};

#endif