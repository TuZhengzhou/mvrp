#include <algorithm>
#include <assert.h>
#include <map>
#include <iostream>
#include <iterator>
#include "libff/common/utils.hpp"
#include "libff/algebra/curves/public_params.hpp"
#include "structs.hpp"
#include "circuit_prover.hpp"
#include "linear_combination.hpp"

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

ProverCircuit::ProverCircuit(const CCredSRS& srs, const vector<LinearCombination>& constraints, \
                const size_t num_multi_gate, const size_t num_commit_input, const map<size_t, size_t>& map, \
                const CSet& S, const vector<Fr>& aL, const vector<Fr>& aR, const vector<Fr>& aO, \
                const IPAProveSystem& ipa_sys) 
{
  // 公开输入的复制
  _srs = CCredSRS(srs);
  _constraints = constraints;
  _map = map;

  // 秘密输入的复制
  _S = CSet(S);
  _a_L = aL;
  _a_R = aR;
  _a_O = aO;

  // 公开输入的推断
  _N = _srs.N;
  _num_constraints = constraints.size();
  _num_multi_gate = num_multi_gate;
  _num_commit_input = num_commit_input;

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
  for(size_t i = 1; i <= _num_multi_gate; i++) {
    _beta_powers_sum_neg = _beta_powers_sum_neg - _srs.get_g1_beta_exp((_N+1)+i);
  }
}

/// Use a challenge, `z`, to flatten the constraints in the
/// constraint system into vectors used for proving and
/// verification.
///
/// # Output
///
/// Returns a tuple of
/// ```text
/// (wL, wR, wO)
/// ```
/// where `w{L,R,O}` is \\( z \cdot z^Q \cdot W_{L,R,O} \\).
void ProverCircuit::flatten_constraints(const vector<LinearCombination>& constraints, const Fr& z, size_t num_multi_gate, vector<Fr>& wL, vector<Fr>& wR, vector<Fr>& wO) {
  
  wL.resize(num_multi_gate, Fr(0));
  wR.resize(num_multi_gate, Fr(0));
  wO.resize(num_multi_gate, Fr(0));

  Fr exp_z = z;
  for(auto lc: constraints) {
    for(auto term: lc.terms) {
      switch (term.first.type)
      {
      case VType::MultiplierLeft:
        wL[term.first.index] += term.second * exp_z;
        break;

      case VType::MultiplierRight:
        wR[term.first.index] += term.second * exp_z;
        break;

      case VType::MultiplierOutput:
        wO[term.first.index] += term.second * exp_z;
        break;

      case VType::Committed:
        break;
      case VType::One:
        break;
      default:
        break;
      }
    }
    exp_z = exp_z * z;
  }

  return;
}

CCircuitProof ProverCircuit::prove(Agenda& agenda, const bool improved) {
  // if(improved) 
  //   return prove_improved(agenda);
  return prove_base(agenda);
}

CCircuitProof ProverCircuit::prove_base(Agenda& agenda) {

  agenda.create_item("Prove_base");
  
  CCircuitProof pi; // global proof pi
  std::ostream_iterator<Fr>   outFr(cout, ", ");
  std::ostream_iterator<bool> outBo(cout, ", ");

  vector<Fr> transcript_fr;
  vector<G1> transcript_g1;
  vector<G2> transcript_g2;
  vector<GT> transcript_gt;
  vector<Fr> random_challenges;

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
    _K_L = std::vector<Fr>(_num_multi_gate);
    _K_R = std::vector<Fr>(_num_multi_gate);
    // Fr::random_e
    for(i = 0; i < _num_multi_gate; i++) {
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
    for(size_t i = 1; i <= _num_multi_gate; i++) {
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
  std::vector<Fr> flatten_wL(_num_multi_gate, Fr(0));
  std::vector<Fr> flatten_wR(_num_multi_gate, Fr(0));
  std::vector<Fr> flatten_wO(_num_multi_gate, Fr(0));

  // G1 point_proof;   // point proof of one item
  transcript_g1.push_back(pi._AI);
  transcript_g1.push_back(pi._AO);
  transcript_g1.push_back(pi._K);
  random_challenges = gen_random_field_elements_from_transcripts(transcript_g1, transcript_g2, transcript_gt, 2);
  _y = random_challenges[0];
  _z = random_challenges[1];

  z_powers_0 = vector_powers(_z, _num_constraints); // 1, z, z^2, ... , z^m, z^(m+1)
  z_powers_1.resize(_num_constraints);
  for(size_t i = 0; i < _num_constraints; i++) {
    z_powers_1[i] = z_powers_0[i+1];
  }
  y_powers_0 = vector_powers(_y, _num_multi_gate - 1);
  y_inv_powers_0 = vector_powers(_y.inverse(), _num_multi_gate - 1);

  flatten_constraints(_constraints, _z, _num_multi_gate, flatten_wL, flatten_wR, flatten_wO);

  vector<Fr> l_x_3;
  vector<Fr> r_x_1, r_x_2, r_x_3, r_x_4, r_x_5;
  l_x_3 = hadmard_product(y_inv_powers_0, flatten_wR);
  r_x_3 = flatten_wL;
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
    r_x_4 = flatten_wO;
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

  for(size_t i = 1; i <= this->_num_multi_gate; i++) {
    g_vec.push_back(this->_srs.get_g1_beta_exp(i));
    h_vec.push_back(y_inv_powers_0[i-1] * this->_srs.get_g1_beta_exp((this->_N+1)+i));
  }

  vector<Fr> bold_wL, bold_wR, bold_wO;
  G1 WL_apos, WR_apos, WO_apos;
  G1 P_G1, P_G1_apos;
  G1 P_G1_v, P_G1_apos_v;

  const bool round_3 = true;
  if(round_3) {

#ifdef CRED_DEBUG
    libff::enter_block("/** Round3: range_prove, output \\mu, r_x, pi_ipa, t_tilde = t(x) **/");
#endif

    transcript_g1.push_back(pi._pi_tilde);
    transcript_gt.push_back(pi._T_1);
    transcript_gt.push_back(pi._T_3);
    transcript_gt.push_back(pi._T_4);
    transcript_gt.push_back(pi._T_5);
    transcript_gt.push_back(pi._T_6);
    _x = gen_random_field_elements_from_transcripts(transcript_g1, transcript_g2, transcript_gt, 1)[0];

    _l_bold = vector<Fr>(_num_multi_gate, Fr(0));
    _r_bold = vector<Fr>(_num_multi_gate, Fr(0));
    for(auto l_coef: l_coefs) {
      _l_bold = vector_add(_l_bold, numerical_mult(_x^(l_coef.second), l_coef.first));
    }
    for(auto r_coef: r_coefs) {
      _r_bold = vector_add(_r_bold, numerical_mult(_x^(r_coef.second), r_coef.first));
    }
    _t_tilde = Fr(0);
    for(auto t_coef: map_t_coefs_6){
      _t_tilde += (_x ^ (t_coef.first)) *  t_coef.second;
    }

    _theta_x = Fr(0);
    for(size_t i = 1; i <= 6; i++) {
      if(i != 2){
        _theta_x += map_theta_6[i] * (_x ^ i);
      }
    }
    _mu = _r_1 * _x  + _r_2 * (_x^2) + _r_3 * (_x^3);

    pi._theta_x = _theta_x;
    pi._mu      = _mu;
    pi._t_tilde = _t_tilde;

    /** round 4 */
    transcript_fr.push_back(pi._theta_x);
    transcript_fr.push_back(pi._mu);
    transcript_fr.push_back(pi._t_tilde);
    _x_t_tilde = gen_random_field_elements_from_transcripts(transcript_fr, transcript_g1, transcript_g2, transcript_gt, (size_t)1)[0];

    libff::enter_block("compute wL, wR, wO");
    bold_wL = hadmard_product(y_inv_powers_0, flatten_wL);
    bold_wR = hadmard_product(y_inv_powers_0, flatten_wR);
    bold_wO = hadmard_product(y_inv_powers_0, flatten_wO);

    WL_apos = G1::zero();
    WR_apos = G1::zero();
    WO_apos = G1::zero();
    for(size_t i = 1; i <= _num_multi_gate; i++) {
      WL_apos = WL_apos + bold_wL[i-1] * _srs.get_g1_beta_exp((_N+1)+i);
      WR_apos = WR_apos + bold_wR[i-1] * _srs.get_g1_beta_exp(i);
      WO_apos = WO_apos + bold_wO[i-1] * _srs.get_g1_beta_exp((_N+1)+i);
    }
    libff::leave_block("compute wL, wR, wO");

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
    for(size_t i = 0; i < _num_multi_gate; i++) {
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
    libff::enter_block("Prove::IpaProve");
    auto pi_ipa             = this->_ipa_prove_sys.IpaProve(g_vec, h_vec, _l_bold, _r_bold, _t_tilde);
    libff::leave_block("Prove::IpaProve");
    agenda.mark_item_end("Prove::IpaProve");

#ifdef CRED_DEBUG
    libff::leave_block("/** Round3: range_prove, output \\mu, r_x, pi_ipa, t_tilde = t(x) **/");
#endif

    pi._pi_ipa = pi_ipa;
  }

  agenda.mark_item_end("Prove_base::Round3");

  vector<cred::Fr> poly_fwr_apos;
  poly_fwr_apos = shift_and_reverse(bold_wR, 1);
  pi._pi_kzg_fwr_apos = kzg_prove(_srs, poly_fwr_apos, poly_fwr_apos.size(), pi._commit_fwr_apos);

  vector<cred::Fr> poly_fwlo_apos;
  poly_fwlo_apos = shift_and_reverse(vector_add(numerical_mult(_x, bold_wL), bold_wO), _N+2);
  pi._pi_kzg_fwlo_apos = kzg_prove(_srs, poly_fwlo_apos, poly_fwlo_apos.size(), pi._commit_fwlo_apos);

  /** compute randomes used for inner-product-argument folding */
  size_t recursion_time;
  std::vector<Fr> randoms; 

  recursion_time = pi._pi_ipa.L_vec.size();
  randoms = std::vector<Fr>(recursion_time, Fr::zero());

  randoms[0] = IPAProveSystem::generate_random(Fr::zero(), pi._pi_ipa.L_vec[0], pi._pi_ipa.R_vec[0]);
  for(size_t i = 1; i < recursion_time; i++) {
    randoms[i] = IPAProveSystem::generate_random(randoms[i-1], pi._pi_ipa.L_vec[i], pi._pi_ipa.R_vec[i]);
  }

  /** compute polynomial fg, fh */
  std::vector<Fr> ss, ss_inv;

  agenda.create_item("Prove::multi_exponentiation_n");
  ss = multi_exponentiation_n(randoms, recursion_time);
  ss_inv.resize(ss.size());

  for(i = 0; i < ss.size(); i++) {
    ss_inv[i] = ss[i].inverse();
  }
  agenda.mark_item_end("Prove::multi_exponentiation_n");

  std::vector<Fr> poly_fg, poly_fh;
  poly_fg = shift_and_reverse(ss, 1);
  pi._pi_kzg_fg = kzg_prove(_srs, poly_fg, poly_fg.size(), pi._commit_fg);

  poly_fh = shift_and_reverse(hadmard_product(y_inv_powers_0, ss_inv), _N+2);
  pi._pi_kzg_fh = kzg_prove(_srs, poly_fh, poly_fh.size(), pi._commit_fh);

  agenda.mark_item_end("Prove_base");
  agenda.mark_mem_usage("Prove_base");

  return pi;
};

} // namespace cred