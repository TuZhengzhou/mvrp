#include <algorithm>
#include <assert.h>
#include <map>
#include <iostream>
#include <iterator>
#include "libff/common/utils.hpp"
#include "libff/algebra/curves/public_params.hpp"
#include "structs.hpp"
#include "verifier_circuit.hpp"
#include "kzg.hpp"

namespace cred {
VerifierCircuit::VerifierCircuit(const CredSRS& srs, const G1& V, \
  const vector<vector<Fr>>& WL, const vector<vector<Fr>>& WR, const vector<vector<Fr>>& WO, const vector<vector<Fr>>& WU, \
  vector<Fr> c, map<size_t, size_t> map, \
  const IPAProveSystem& ipa_sys) 
{
// 公开输入的复制
  _srs = CredSRS(srs);
  _V  = G1(V);
  _WL = WL;
  _WR = WR;
  _WO = WO;
  _WU = WU;
  _c = c;
  _map = map;

  // 公开输入的推断
  _N = _srs.N;
  _Q = WL.size();
  _n = WL[0].size();
  _num_commit_input = WU[0].size();

  // 子证明系统
  _ipa_prove_sys = ipa_sys;

  _beta_powers_sum_neg = G1::zero();
  for(size_t i = 1; i <= _n; i++) {
    _beta_powers_sum_neg = _beta_powers_sum_neg - _srs.get_g1_beta_exp((_N+1)+i);
  }
}

bool VerifierCircuit::verify(const ProofCircuit& pi, Agenda& agenda, const bool improved) {
  // if(improved)
    // return verify_improved(pi, agenda);
  return verify_base(pi, agenda);
}

bool VerifierCircuit::verify_base(const ProofCircuit& pi, Agenda& agenda) {

  bool result;
  agenda.create_item("Verify_base");

  vector<Fr> vec_fr;
  vector<G1> vec_1;
  vector<G2> vec_2;
  vector<GT> vec_T;
  vector<Fr> vec_r;

  GT beta_N_plus_2_GT = ReducedPairing(this->_srs.get_g1_beta_exp(_N), this->_srs.get_g2_beta_exp(2));
  Fr y, z;      // randomes for second round
  Fr x;         // randome for third round
  Fr x_t_tilde; // random for last round

  Fr delta;
  std::vector<Fr> z_powers_0; // index 0 with z^0, [1, z, z^2, ... , z^m, z^(m+1)]
  std::vector<Fr> z_powers_1; // index 0 with z^0, [1, z, z^2, ... , z^m, z^(m+1)]
  std::vector<Fr> y_powers_0; // index 0 with y^1, [1, y, y^2, ..., y^(D-1)]
  std::vector<Fr> y_inv_powers_0;

  // select random y and z
  vec_1.push_back(pi._AI);
  vec_1.push_back(pi._AO);
  vec_1.push_back(pi._K);
  vec_r = generate_random_fr_vec(vec_1, vec_2, vec_T, 2);
  y = vec_r[0];
  z = vec_r[1];

  // precompute some data
  z_powers_0 = vector_powers(z, _Q); // 1, z, z^2, ... , z^m, z^(m+1)
  z_powers_1.resize(_Q);
  for(size_t i = 0; i < _Q; i++) {
    z_powers_1[i] = z_powers_0[i+1];
  }
  y_powers_0 = vector_powers(y, _n - 1);
  y_inv_powers_0 = vector_powers(y.inverse(), _n - 1);

  vector<Fr> l_x_3;
  vector<Fr> r_x_3;
  l_x_3 = hadmard_product(y_inv_powers_0, vec_matrix_mult(z_powers_1, _WR, true));
  r_x_3 = vec_matrix_mult(z_powers_1, _WL, true);
  delta = inner_product(l_x_3, r_x_3);

  // select random x
  vec_1.push_back(pi._pi_tilde);
  vec_T.push_back(pi._T_1);
  vec_T.push_back(pi._T_3);
  vec_T.push_back(pi._T_4);
  vec_T.push_back(pi._T_5);
  vec_T.push_back(pi._T_6);
  x = generate_random_fr_vec(vec_1, vec_2, vec_T, 1)[0];

  // select random x_t_tilde
  vec_fr.push_back(pi._theta_x);
  vec_fr.push_back(pi._mu);
  vec_fr.push_back(pi._t_tilde);
  x_t_tilde = generate_random_fr_vec(vec_fr, vec_1, vec_2, vec_T, (size_t)1)[0];

  // another pre computations
  vector<Fr> bold_wL, bold_wR, bold_wO;
  G1 WL_apos, WR_apos, WO_apos;
  G1 P_G1, P_G1_apos;
  G1 P_G1_v, P_G1_apos_v;

  bold_wL = hadmard_product(y_inv_powers_0, vec_matrix_mult(z_powers_1, _WL, true));
  bold_wR = hadmard_product(y_inv_powers_0, vec_matrix_mult(z_powers_1, _WR, true));
  bold_wO = hadmard_product(y_inv_powers_0, vec_matrix_mult(z_powers_1, _WO, true));

  WL_apos = G1::zero();
  WR_apos = G1::zero();
  WO_apos = G1::zero();
  for(size_t i = 1; i <= _n; i++) {
    WL_apos = WL_apos + bold_wL[i-1] * _srs.get_g1_beta_exp((_N+1)+i);
    WR_apos = WR_apos + bold_wR[i-1] * _srs.get_g1_beta_exp(i);
    WO_apos = WO_apos + bold_wO[i-1] * _srs.get_g1_beta_exp((_N+1)+i);
  }

  P_G1_v = G1::zero();
  P_G1_v = P_G1_v + (x ^ 1) * (pi._AI);
  P_G1_v = P_G1_v + (x ^ 1) * (WL_apos + WR_apos);
  P_G1_v = P_G1_v + (x ^ 2) * pi._AO;
  P_G1_v = P_G1_v + (x ^ 3) * pi._K;
  P_G1_v = P_G1_v + WO_apos;
  P_G1_v = P_G1_v + (-pi._mu) * _srs.get_g1_beta_exp(_N);
  P_G1_v = P_G1_v + _beta_powers_sum_neg;
  
  // this->_ipa_prove_sys._u = x_t_tilde * _srs.get_g1_beta_exp(_N);
  // P_G1_apos = P_G1 + pi._t_tilde * this->_ipa_prove_sys._u;

#ifdef CRED_DEBUG
  libff::enter_block("/** varient of pointproof verify */");
#endif

  // first check
  /** compute
   * [left_1]_T = (x · [T1]_T )
   * [left_2]_T = (x^2 · [T1]_T )
   * [left_3]_T = e( C, [∑j in I {z^(1+j) β^(N+1−j)}]_2 )
   * [left]_T = [left_1]_T * [left_2]_T * [left_3]_T
  **/
  G2 sum_z_pows_beta_pows; // sum_z_pows_beta_pows = [∑j in I {z^(1+j) β^(N+1−j)}]_2

  sum_z_pows_beta_pows = G2::zero();
  for(size_t j = 1; j <= _num_commit_input; j++) {
    sum_z_pows_beta_pows = sum_z_pows_beta_pows + z_powers_0[j] * this->_srs.get_g2_beta_exp(_N + 1 - _map[j]);
  }

  GT check_l_1, check_l_2;
  GT check_r_1, check_r_2, check_r_3;
  GT check_l, check_r;

  check_l_1 = (pi._T_1) ^ x;
  check_l_1 = check_l_1 * ((pi._T_3) ^ (x^3));
  check_l_1 = check_l_1 * ((pi._T_4) ^ (x^4));
  check_l_1 = check_l_1 * ((pi._T_5) ^ (x^5));
  check_l_1 = check_l_1 * ((pi._T_6) ^ (x^6));
  check_l_2 = ReducedPairing(_V, (x^2) * sum_z_pows_beta_pows);

  check_r_1 = _srs.gt ^ (pi._t_tilde - (x^2) * (delta + inner_product(z_powers_1, _c)));
  check_r_2 = beta_N_plus_2_GT ^ (pi._theta_x);
  check_r_3 = ReducedPairing((x^2) * pi._pi_tilde, _srs.g2_base_);

  check_l = check_l_1 * check_l_2;
  check_r = check_r_1 * check_r_2 * check_r_3;

  result = (check_l == check_r);

#ifdef CRED_DEBUG
  libff::leave_block("/** varient of pointproof verify */");
  assert(result);
  libff::enter_block("/** inner product check for range proof **/");
#endif

  agenda.create_item("/** inner product check for range proof **/");

  std::vector<G1> g_vec;          // [[β^1]_1, ..., [β^D]_1]
  std::vector<G1> h_vec;          // [[y^(1-1) * β^(N+1+1)]_1, ..., [[y^(1-i) * β^(N+1+i)]_1], ..., [y^(1-D) * β^(N+1+D)]_1]

  for(size_t i = 1; i <= this->_n; i++) {
    g_vec.push_back(this->_srs.get_g1_beta_exp(i));
    h_vec.push_back(y_inv_powers_0[i-1] * this->_srs.get_g1_beta_exp((this->_N+1)+i));
  }

  this->_ipa_prove_sys._u = x_t_tilde * _srs.get_g1_beta_exp(_N);
  this->_ipa_prove_sys._P = P_G1_v;
  this->_ipa_prove_sys._c = pi._t_tilde;
  // P_G1_apos = P_G1 + pi._t_tilde * this->_ipa_prove_sys._u;

  agenda.create_item("Verify_base::IpaVerify");
  result &= this->_ipa_prove_sys.IpaMultiExpVerify(pi._pi_ipa, g_vec, h_vec);
  agenda.mark_item_end("Verify_base::IpaVerify");
  agenda.mark_item_end("/** inner product check for range proof **/");

  
#ifdef CRED_DEBUG
    libff::leave_block("/** inner product check for range proof **/");
#endif

#ifdef CRED_DEBUG
  assert(result);
#endif

  agenda.mark_item_end("Verify_base");
  agenda.mark_mem_usage("Verify_base");

  return result;
}

}