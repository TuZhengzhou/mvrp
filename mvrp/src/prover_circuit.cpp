#include <algorithm>
#include <assert.h>
#include <map>
#include <iostream>
#include <iterator>
#include "libff/common/utils.hpp"
#include "libff/algebra/curves/public_params.hpp"
#include "structs.hpp"
#include "prover_circuit.hpp"
#include "kzg.hpp"

namespace cred {

void ProverCircuit::commit_set() {
  this->_r = Fr::random_element();
  _V = this->_r * this->_srs.get_g1_beta_exp(_N);  // item with index n havs exponent n+1
  for(size_t i = 1; i <= _S.m; i++) {
    _V = _V + _S.get_set_value(i) * this->_srs.get_g1_beta_exp(i);
  }
}

void ProverCircuit::point_proof_pre_compute() {

  _point_proofs.resize(_S.m);
  
  size_t i, j;
  G1 point_proof;

  for(j = 1; j <= _S.m; j++) {
    point_proof = this->_r * _srs.get_g1_beta_exp(2 * _N + 1 - j);
    for(i = 1; i <= _S.m; i++) {
      if(i != j) {
        point_proof = point_proof + _S.get_set_value(i) * this->_srs.get_g1_beta_exp(_N + 1 - j + i);
      }
    }
    _point_proofs[j-1] = point_proof;
  }
}

const G1& ProverCircuit::point_proof(size_t idx) {
#ifdef CRED_DEBUG
  assert(idx >= 1 && idx <= _S.set_values.size());
#endif
  return _point_proofs[idx-1];
}

// ProverCircuit::ProverCircuit(const CredSRS& srs, const SET& S, const Ranges& ranges, const IPAProveSystem& ipa_sys) {
//   _srs = CredSRS(srs);
//   _S = SET(S);
//   _ipa_prove_sys = ipa_sys;

//   _N = _srs.N;

//   commit_set();                 // 
//   point_proof_pre_compute();
// }
ProverCircuit::ProverCircuit(const CredSRS& srs, \
                const vector<vector<Fr>>& WL, const vector<vector<Fr>>& WR, const vector<vector<Fr>>& WO, const vector<vector<Fr>>& WU, \
                vector<Fr> c, map<size_t, size_t> map, \
                const SET& S, vector<Fr> aL, vector<Fr> aR, vector<Fr> aO, \
                const IPAProveSystem& ipa_sys) 
{
  // 公开输入的复制
  _srs = CredSRS(srs);
  _WL = WL;
  _WR = WR;
  _WO = WO;
  _WU = WU;
  _c = c;
  _map = map;

  // 秘密输入的复制
  _S = SET(S);
  _a_L = aL;
  _a_R = aR;
  _a_O = aO;

  // 公开输入的推断
  _N = _srs.N;
  _Q = WL.size();
  _n = WL[0].size();
  _num_commit_input = WU[0].size();

  // 秘密输入的推断
  _u.clear();
  for(size_t i = 1; i <= _num_commit_input; i++) {
    _u.push_back(_S.get_set_value(_map[i]));
  }

  // 子证明系统
  _ipa_prove_sys = ipa_sys;

  commit_set();                 // 
  point_proof_pre_compute();

  _beta_powers_sum_neg = G1::zero();
  for(size_t i = 1; i <= _n; i++) {
    _beta_powers_sum_neg = _beta_powers_sum_neg - _srs.get_g1_beta_exp((_N+1)+i);
  }
}

ProofCircuit ProverCircuit::prove(Agenda& agenda, const bool improved) {
  // if(improved) 
  //   return prove_improved(agenda);
  return prove_base(agenda);
}

ProofCircuit ProverCircuit::prove_base(Agenda& agenda) {

  agenda.create_item("Prove_base");
  
  ProofCircuit pi; // global proof pi
  std::ostream_iterator<Fr>   outFr(cout, ", ");
  std::ostream_iterator<bool> outBo(cout, ", ");

  vector<Fr> vec_fr;
  vector<G1> vec_1;
  vector<G2> vec_2;
  vector<GT> vec_T;
  vector<Fr> vec_r;

  GT beta_N_plus_2_GT = ReducedPairing(this->_srs.get_g1_beta_exp(_N), this->_srs.get_g2_beta_exp(2));

  // round 1
  Fr _r_1, _r_2, _r_3;
  std::vector<Fr> _K_L, _K_R;
  G1 _AI, _AO, _K;

  // round 2
  Fr _y, _z;
  Fr delta;
  vector<pair<vector<Fr>, size_t>> l_coefs;
  vector<pair<vector<Fr>, size_t>> r_coefs;
  map<size_t, Fr> map_t_coefs_6;
  map<size_t, Fr> map_theta_6;
  map<size_t, GT> map_T_6;
  G1 _pi_tilde;

  // round 3
  Fr _x, _t_tilde, _theta_x, _mu;
  std::vector<Fr> _l_bold, _r_bold;

  // round 4
  Fr _x_t_tilde;
  
  agenda.create_item("Prove_base::Round1");
  /**
    Round1: commit_for_ipa, output [AI]_1, [AO]_1 [K]_1
  **/
  size_t base, i;

  const bool round_1 = true;
  if(round_1) {
    /** sample r1, r2 ←$ Zp, KL, KR ←$ Z^D_p **/
    _r_1 = Fr::random_element();
    _r_2 = Fr::random_element();
    _r_3 = Fr::random_element();
    _K_L = std::vector<Fr>(_n);
    _K_R = std::vector<Fr>(_n);
    // Fr::random_e
    for(i = 0; i < _n; i++) {
      _K_L[i] = Fr::random_element();
      _K_R[i] = Fr::random_element();
    }

    /** compute 
     * [AI]1 = [D∑i=1(aL[i]β^i + aR[i]β^(N+1+i)) + r_1 β^N ]_1 and 
     * [AO]1 = [D∑i=1(aO[i]β^i) + r_2 β^N ]_1 and 
     * [K]1 = [D∑i=1(KL[i]β^i + KR[i]β^(N+1+i)) + r_3 β^N ]_1
    **/
#ifdef CRED_DEBUG
    libff::enter_block("/** compute [A]_1, [K]_1 **/");
#endif

    _AI = _r_1 * this->_srs.get_g1_beta_exp(_N);
    _AO = _r_2 * this->_srs.get_g1_beta_exp(_N);
    _K  = _r_3 * this->_srs.get_g1_beta_exp(_N);
    for(size_t i = 1; i <= _n; i++) {
      _AI = _AI + _a_L[i-1] * this->_srs.get_g1_beta_exp(i) + _a_R[i-1] * this->_srs.get_g1_beta_exp((_N+1)+i);
      _AO = _AO + _a_O[i-1] * this->_srs.get_g1_beta_exp(i);
      _K = _K + _K_L[i-1] * this->_srs.get_g1_beta_exp(i) + _K_R[i-1] * this->_srs.get_g1_beta_exp((_N+1)+i);
    }
#ifdef CRED_DEBUG
    libff::leave_block("/** compute [A]_1, [K]_1 **/");
#endif
    pi._AI = _AI;
    pi._AO = _AO;
    pi._K = _K;
  }
  agenda.mark_item_end("Prove_base::Round1");

  agenda.create_item("Prove_base::Round2");
  /**
    Round2: point_prove, output [pi_tilde]_1, [T_1]_T, [T_2]_T
  **/
  std::vector<Fr> z_powers_0; // index 0 with z^0, [1, z, z^2, ... , z^m, z^(m+1)]
  std::vector<Fr> z_powers_1; // index 0 with z^0, [1, z, z^2, ... , z^m, z^(m+1)]
  std::vector<Fr> y_powers_0; // index 0 with y^1, [1, y, y^2, ..., y^(D-1)]
  std::vector<Fr> y_inv_powers_0;

  // G1 point_proof;   // point proof of one item
  vec_1.push_back(pi._AI);
  vec_1.push_back(pi._AO);
  vec_1.push_back(pi._K);
  vec_r = generate_random_fr_vec(vec_1, vec_2, vec_T, 2);
  _y = vec_r[0];
  _z = vec_r[1];

  z_powers_0 = vector_powers(_z, _Q); // 1, z, z^2, ... , z^m, z^(m+1)
  z_powers_1.resize(_Q);
  for(size_t i = 0; i < _Q; i++) {
    z_powers_1[i] = z_powers_0[i+1];
  }
  y_powers_0 = vector_powers(_y, _n - 1);
  y_inv_powers_0 = vector_powers(_y.inverse(), _n - 1);

  vector<Fr> l_x_3;
  vector<Fr> r_x_1, r_x_2, r_x_3, r_x_4, r_x_5;
  l_x_3 = hadmard_product(y_inv_powers_0, vec_matrix_mult(z_powers_1, _WR, true));
  r_x_3 = vec_matrix_mult(z_powers_1, _WL, true);
  delta = inner_product(l_x_3, r_x_3);
  
  const bool round_2 = true;
  if(round_2) {

#ifdef CRED_DEBUG
    libff::enter_block("/** Round2: point_prove, output [pi_tilde]_1, [T_1]_T, [T_2]_T **/");
#endif

    l_coefs.emplace_back(pair<vector<Fr>, size_t>(_a_L , (size_t)1));
    l_coefs.emplace_back(pair<vector<Fr>, size_t>(_a_O , (size_t)2));
    l_coefs.emplace_back(pair<vector<Fr>, size_t>(l_x_3, (size_t)1));
    l_coefs.emplace_back(pair<vector<Fr>, size_t>(_K_L , (size_t)3));
    r_x_1 = hadmard_product(y_powers_0, _a_R);
    r_x_2 = vector_neg(y_powers_0);
    r_x_4 = vec_matrix_mult(z_powers_1, _WO, true);
    r_x_5 = hadmard_product(y_powers_0, _K_R);
    r_coefs.emplace_back(pair<vector<Fr>, size_t>(r_x_1 , (size_t)1));
    r_coefs.emplace_back(pair<vector<Fr>, size_t>(r_x_2 , (size_t)0));
    r_coefs.emplace_back(pair<vector<Fr>, size_t>(r_x_3 , (size_t)1));
    r_coefs.emplace_back(pair<vector<Fr>, size_t>(r_x_4 , (size_t)0));
    r_coefs.emplace_back(pair<vector<Fr>, size_t>(r_x_5 , (size_t)3));

    for(auto l_coef: l_coefs) {
      for(auto r_coef: r_coefs) {
        auto key = l_coef.second + r_coef.second;
        if(map_t_coefs_6.find(key) == map_t_coefs_6.end()){
          map_t_coefs_6[key] = inner_product(l_coef.first, r_coef.first);
        } else {
          map_t_coefs_6[key] += inner_product(l_coef.first, r_coef.first);
        }
      }
    }
    for(size_t i = 1; i <= 6; i++) {
      map_theta_6[i] = Fr::random_element();
      map_T_6[i] = (_srs.gt ^ map_t_coefs_6[i]) * (beta_N_plus_2_GT ^ map_theta_6[i]);
    }

    _pi_tilde = G1::zero();
    for(size_t j = 1; j <= _num_commit_input; j++) {
      _pi_tilde = _pi_tilde + z_powers_0[j] * this->point_proof(_map[j]);
    }

    pi._T_1 = map_T_6[1];
    pi._T_3 = map_T_6[3];
    pi._T_4 = map_T_6[4];
    pi._T_5 = map_T_6[5];
    pi._T_6 = map_T_6[6];
    pi._pi_tilde = _pi_tilde;
  }
  agenda.mark_item_end("Prove_base::Round2");
#ifdef CRED_DEBUG
    libff::leave_block("/** Round2: point_prove, output [pi_tilde]_1, [T_1]_T, [T_2]_T **/");
#endif

  agenda.create_item("Prove_base::Round3");
  /**
    Round3: range_prove, output \mu, r_x, pi_ipa, t_tilde = t(x)
  **/
  std::vector<G1> g_vec;          // [[β^1]_1, ..., [β^D]_1]
  std::vector<G1> h_vec;          // [[y^(1-1) * β^(N+1+1)]_1, ..., [[y^(1-i) * β^(N+1+i)]_1], ..., [y^(1-D) * β^(N+1+D)]_1]

  for(size_t i = 1; i <= this->_n; i++) {
    g_vec.push_back(this->_srs.get_g1_beta_exp(i));
    h_vec.push_back(y_inv_powers_0[i-1] * this->_srs.get_g1_beta_exp((this->_N+1)+i));
  }

  const bool round_3 = true;
  if(round_3) {

#ifdef CRED_DEBUG
    libff::enter_block("/** Round3: range_prove, output \\mu, r_x, pi_ipa, t_tilde = t(x) **/");
#endif

    vec_1.push_back(pi._pi_tilde);
    vec_T.push_back(pi._T_1);
    vec_T.push_back(pi._T_3);
    vec_T.push_back(pi._T_4);
    vec_T.push_back(pi._T_5);
    vec_T.push_back(pi._T_6);
    _x = generate_random_fr_vec(vec_1, vec_2, vec_T, 1)[0];

    _l_bold = vector<Fr>(_n, Fr::zero());
    _r_bold = vector<Fr>(_n, Fr::zero());
    for(auto l_coef: l_coefs) {
      _l_bold = vector_add(_l_bold, numerical_mult(_x^(l_coef.second), l_coef.first));
    }
    for(auto r_coef: r_coefs) {
      _r_bold = vector_add(_r_bold, numerical_mult(_x^(r_coef.second), r_coef.first));
    }
    _t_tilde = Fr::zero();
    for(auto t_coef: map_t_coefs_6){
      _t_tilde += (_x ^ (t_coef.first)) *  t_coef.second;
    }

    _theta_x = Fr::zero();
    for(size_t i = 1; i <= 6; i++) {
      if(i != 2){
        _theta_x += map_theta_6[i] * (_x ^ i);
      }
    }
    _mu = _r_1 * _x  + _r_2 * (_x^2) + _r_3 * (_x^3);

    pi._theta_x = _theta_x;
    pi._mu      = _mu;
    pi._t_tilde = _t_tilde;

#ifdef CRED_DEBUG
    assert(_t_tilde == inner_product(_l_bold, _r_bold));

    Fr theta_sum, t_sum;
    GT check_l_1, check_t;

    theta_sum = map_theta_6[1] * _x;
    check_t   = _srs.gt ^ (map_t_coefs_6[1] * _x);
    check_l_1 = map_T_6[1] ^ _x;
    for(size_t i = 2; i <= 6; i++) {
      theta_sum = theta_sum + map_theta_6[i] * (_x^i);
      check_t   = check_t   * (_srs.gt ^ (map_t_coefs_6[i] * (_x ^ i)));
      check_l_1 = check_l_1 * ((map_T_6[i]) ^ (_x^i));
    }
    assert((_srs.gt ^ (_t_tilde)) == check_t);
    assert((_srs.gt ^ (_t_tilde)) * (beta_N_plus_2_GT ^ theta_sum) == check_l_1);

    vector<Fr> w_tmp = vector<Fr>(_Q, Fr::zero());
    w_tmp = vector_add(w_tmp, vec_matrix_mult(_a_L, _WL, false));
    w_tmp = vector_add(w_tmp, vec_matrix_mult(_a_R, _WR, false));
    w_tmp = vector_add(w_tmp, vec_matrix_mult(_a_O, _WO, false));
    assert(map_t_coefs_6[2] == (inner_product(z_powers_1, w_tmp) + delta));  // delta 计算正确

    G2 sum_z_pows_beta_pows = G2::zero();
    for(size_t j = 1; j <= _num_commit_input; j++) {
      sum_z_pows_beta_pows = sum_z_pows_beta_pows + z_powers_0[j] * this->_srs.get_g2_beta_exp(_N + 1 - _map[j]);
    }

    GT check_l_2;
    GT check_r_1, check_r_2, check_r_3;
    GT check_l, check_r;

    check_l_1 = (pi._T_1) ^ _x;
    check_l_1 = check_l_1 * ((pi._T_3) ^ (_x^3));
    check_l_1 = check_l_1 * ((pi._T_4) ^ (_x^4));
    check_l_1 = check_l_1 * ((pi._T_5) ^ (_x^5));
    check_l_1 = check_l_1 * ((pi._T_6) ^ (_x^6));
    check_l_2 = ReducedPairing(_V, (_x^2) * sum_z_pows_beta_pows);

    assert(inner_product(z_powers_1, _c) == Fr::zero());  // c 计算正确

    check_r_1 = _srs.gt ^ (pi._t_tilde - (_x^2) * (delta + inner_product(z_powers_1, _c)));
    check_r_2 = beta_N_plus_2_GT ^ (pi._theta_x);
    check_r_3 = ReducedPairing((_x^2) * pi._pi_tilde, _srs.g2_base_);

    check_l = check_l_1 * check_l_2;
    check_r = check_r_1 * check_r_2 * check_r_3;

    Fr mizi = Fr::zero();
    for (size_t i = 1; i <= _num_commit_input; i++) {
      mizi = mizi + z_powers_0[i] * _S.get_set_value(_map[i]);
    }

    GT mizi_gt = _srs.gt ^ (mizi * (_x^2));

    for(size_t tmp = 1; tmp <= _num_commit_input; tmp++) {
      assert(ReducedPairing(_V, _srs.get_g2_beta_exp(_N+1-_map[tmp])) == ReducedPairing(this->point_proof(_map[tmp]), _srs.g2_base_) * ((_srs.gt)^(_S.get_set_value(_map[tmp]))));
    }
    G2 sum_z_beta = G2::zero();
    for(size_t j = 1; j <= _num_commit_input; j++) {
      sum_z_beta = sum_z_beta + z_powers_0[j] * this->_srs.get_g2_beta_exp(_N + 1 - _map[j]);
    }
    G1 pi_tilde = G1::zero();
    for(size_t j = 1; j <= _num_commit_input; j++) {
      pi_tilde = pi_tilde + z_powers_0[j] * this->point_proof(_map[j]);
    }
    assert(ReducedPairing(_V, sum_z_beta) == ReducedPairing(pi_tilde, _srs.g2_base_) * (_srs.gt ^ mizi));
    

    assert(check_l_2 == mizi_gt * check_r_3);
    assert(check_l_1 * mizi_gt == check_r_1 * check_r_2);

    assert(check_l == check_r);
#endif
    /** round 4 */
    vec_fr.push_back(pi._theta_x);
    vec_fr.push_back(pi._mu);
    vec_fr.push_back(pi._t_tilde);
    _x_t_tilde = generate_random_fr_vec(vec_fr, vec_1, vec_2, vec_T, (size_t)1)[0];

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
    P_G1_v = P_G1_v + (_x ^ 1) * (pi._AI);
    P_G1_v = P_G1_v + (_x ^ 1) * (WL_apos + WR_apos);
    P_G1_v = P_G1_v + (_x ^ 2) * pi._AO;
    P_G1_v = P_G1_v + (_x ^ 3) * pi._K;
    P_G1_v = P_G1_v + WO_apos;
    P_G1_v = P_G1_v + (-pi._mu) * _srs.get_g1_beta_exp(_N);
    P_G1_v = P_G1_v + _beta_powers_sum_neg;
    
#ifdef CRED_DEBUG
    P_G1 = G1::zero();
    for(size_t i = 0; i < _n; i++) {
      P_G1 = P_G1 + _l_bold[i] * g_vec[i];
      P_G1 = P_G1 + _r_bold[i] * h_vec[i];
    }
    assert(P_G1 == P_G1_v);
#endif
    /************************* buller_ipa ***********************/
    this->_ipa_prove_sys._u = _x_t_tilde * _srs.get_g1_beta_exp(_N);
    P_G1_apos = P_G1 + pi._t_tilde * this->_ipa_prove_sys._u;

    this->_ipa_prove_sys._P = P_G1_apos;
    this->_ipa_prove_sys._c = _t_tilde;
    agenda.create_item("Prove::IpaProve");
    auto pi_ipa             = this->_ipa_prove_sys.IpaProve(g_vec, h_vec, _l_bold, _r_bold, _t_tilde);
    agenda.mark_item_end("Prove::IpaProve");

#ifdef CRED_DEBUG
    libff::leave_block("/** Round3: range_prove, output \\mu, r_x, pi_ipa, t_tilde = t(x) **/");
#endif

    pi._pi_ipa = pi_ipa;
  }

  agenda.mark_item_end("Prove_base::Round3");
  agenda.mark_item_end("Prove_base");
  agenda.mark_mem_usage("Prove_base");

  return pi;
};

}