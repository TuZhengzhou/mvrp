#ifndef _CIRCUIT_HPP_
#define _CIRCUIT_HPP_
#include "basic_types.hpp"
#include "ipa.hpp"
#include "kzg.hpp"

namespace cred {

class CCircuitProof {
public:
  // round 1
  G1 _AI, _AO, _K;

  // round 2
  G1 _pi_tilde; // product_(j \in _I) pi_j
  GT _T_1, _T_3, _T_4, _T_5, _T_6;

  // round 3
  Fr _mu;
  Fr _theta_x;
  Fr _t_tilde;  // t(x)
  CIpaProof _pi_ipa;   // 

  G1 _commit_fwr_apos, _commit_fwlo_apos;
  CKzgProof _pi_kzg_fwr_apos, _pi_kzg_fwlo_apos;

  G1 _commit_fg, _commit_fh;
  CKzgProof _pi_kzg_fg;
  CKzgProof _pi_kzg_fh;

  inline const size_t size_base() const {
    return 3 * sizeof(Fr) + 5 * sizeof(GT) + 4 * sizeof(G1) + _pi_ipa.size();
  }

};

};  // namespace cred

#endif // _CIRCUIT_HPP_