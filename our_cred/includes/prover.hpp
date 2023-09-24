#ifndef CRED_PROVER_HPP
#define CRED_PROVER_HPP
#include <vector>
#include <map>
#include <tuple>
#include "libff/algebra/curves/public_params.hpp"
#include "libff/algebra/fields/field_utils.hpp"
#include "libff/common/utils.hpp"
#include "structs.hpp"
#include "ipa.hpp"

namespace cred {

class Prover{
public:
  CredSRS _srs; // srs
  SET _S;   // set S
  G1 _C;    // commitment C
  Fr _r;    // random r for commitment C
  Ranges _ranges;
  IPAProveSystem _ipa_prove_sys;

  size_t _D;    // D = sum_(j in I) n_(j)
  size_t _N;    // _N = srs.N
  map<size_t, size_t> _part_sums; // sum_(i=1)^(j-1) n_i

  // round 1
  Fr _r_1, _r_2;
  std::vector<Fr> _a_L, _a_R, _K_L, _K_R;
  std::map<size_t, std::vector<Fr>> _ds;
  G1 _A, _K;

  // round 2
  Fr t0_1, t0_2, t0_3;
  Fr _y, _z, _t_0, _t_1, _t_2, _r_3, _r_4;
  G1 _pi_tilde;
  GT _T_1, _T_2;

  // round 3
  Fr _x, _t_tilde, _r_x, _mu;
  std::vector<Fr> _l_bold, _r_bold;
  G1 _P_G1;

  Prover() {};
  Prover(const CredSRS srs, const SET S, const Ranges& ranges, const IPAProveSystem& ipa_sys);

  CredProof prove();
  
private:
  const std::vector<Fr>& get_dj(size_t);

  void commit_set();            // select random _r and generate commitment _C
  bool t0_check(const Fr& two, const vector<Fr>& one_vec, const vector<Fr>& z_powers_0, const vector<Fr>& y_powers_1);
  bool easy_check(const std::vector<Fr>& z_exps); // 直接给元素值, 然后检查
};

};
#endif