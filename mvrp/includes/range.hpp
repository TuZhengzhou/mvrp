#ifndef _RANGE_HPP_
#define _RANGE_HPP_
#include "basic_types.hpp"
#include "ipa.hpp"
#include "kzg.hpp"

namespace cred {

class CRangeProof {
public:
  // round 1
  G1 _A, _K;

  // round 2
  G1 _pi_tilde; // product_(j \in _I) pi_j
  GT _T_1, _T_2;
  

  // round 3
  Fr _mu;   // r_1 + r_2 * x
  Fr _r_x;  // r_3 * x + r_4 * x^2
  Fr _t_tilde;  // t(x)
  CIpaProof _pi_ipa;   // 

  G1 _commit_fzy, _commit_fg, _commit_fh;
  CKzgProof _pi_kzg_fzy, _pi_kzg_fg, _pi_kzg_fh;

  inline const size_t size_base() const {
    return 3 * sizeof(Fr) + 2 * sizeof(GT) + 3 * sizeof(G1) + _pi_ipa.size();
  }
  inline const size_t size_improved() const {
    return size_base() + 3 * sizeof(G1) + 3 * 2 * sizeof(G1); // 3 commitments, 3 times (2 G1 in a kzg proof)
  }
};

};  // namespace cred
#endif // _RANGE_HPP_