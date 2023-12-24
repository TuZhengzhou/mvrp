#include <algorithm>
#include <assert.h>
#include <map>
#include <iostream>
#include <iterator>
#include "libff/common/utils.hpp"
#include "libff/algebra/curves/public_params.hpp"
#include "structs.hpp"
#include "circuit_verifier.hpp"
#include "linear_combination.hpp"

namespace cred {
VerifierCircuit::VerifierCircuit(const CCredSRS& srs, const vector<LinearCombination>& constraints, const G1& V, \
  const size_t num_multi_gate, const size_t num_commit_input, \
  const map<size_t, size_t>& map, const IPAProveSystem& ipa_sys) 
{
// 公开输入的复制
  _srs = CCredSRS(srs);
  _constraints = constraints;
  _V  = G1(V);
  _num_constraints = constraints.size();
  _num_multi_gate = num_multi_gate;
  _num_commit_input = num_commit_input;
  _map = map;

  // 公开输入的推断
  _N = _srs.N;


  // 子证明系统
  _ipa_prove_sys = ipa_sys;

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
/// (wL, wR, wO, wc)
/// ```
/// where `w{L,R,O, c}` is \\( z^Q \cdot W_{L,R,O, c} \\).
///
/// This has the same logic as `ProverCS::flattened_constraints()`
/// but also computes the constant terms (which the prover skips
/// because they're not needed to construct the proof).
void VerifierCircuit::flatten_constraints(const vector<LinearCombination>& constraints, const Fr& z, size_t num_multi_gate, Fr& wc, vector<Fr>& wL, vector<Fr>& wR, vector<Fr>& wO) {
  
  wc = Fr(0);
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
        wc = term.second * exp_z;
        break;
      default:
        break;
      }
    }
    exp_z = exp_z * z;
  }

  return;
}


bool VerifierCircuit::verify(const CCircuitProof& pi, Agenda& agenda, const bool improved) {
  // if(improved)
  //   return verify_improved(pi, agenda);
  return verify_base(pi, agenda);
}

bool VerifierCircuit::verify_base(const CCircuitProof& pi, Agenda& agenda) {

  bool result;
  agenda.create_item("Verify_base");

  vector<Fr> transcript_fr;
  vector<G1> transcript_g1;
  vector<G2> transcript_g2;
  vector<GT> transcript_gt;
  vector<Fr> random_challenges;

  GT beta_N_plus_2_GT = ReducedPairing(this->_srs.get_g1_beta_exp(_N), this->_srs.get_g2_beta_exp(2));
  Fr y, z;      // randomes for second round
  Fr x;         // randome for third round
  Fr x_t_tilde; // random for last round

  Fr delta;
  std::vector<Fr> z_powers_0;     // index 0 with z^0, [1, z, z^2, ... , z^m]
  std::vector<Fr> z_powers_1;     // index 0 with z^1, [   z, z^2, ... , z^m]
  std::vector<Fr> y_powers_0;     // index 0 with y^0, [1, y, y^2, ... , y^m]
  std::vector<Fr> y_inv_powers_0; // index 0 with y^0, [1, y^(-1), y^(-2), ... , y^(-m)]
  std::vector<Fr> flatten_wL(_num_multi_gate, Fr(0));
  std::vector<Fr> flatten_wR(_num_multi_gate, Fr(0));
  std::vector<Fr> flatten_wO(_num_multi_gate, Fr(0));
  Fr flatten_wc;

  // select random y and z
  transcript_g1.push_back(pi._AI);
  transcript_g1.push_back(pi._AO);
  transcript_g1.push_back(pi._K);
  random_challenges = gen_random_field_elements_from_transcripts(transcript_g1, transcript_g2, transcript_gt, 2);
  y = random_challenges[0];
  z = random_challenges[1];

  // precompute some data
  z_powers_0 = vector_powers(z, _num_constraints); // 1, z, z^2, ... , z^m, z^(m+1)
  z_powers_1.resize(_num_constraints);
  for(size_t i = 0; i < _num_constraints; i++) {
    z_powers_1[i] = z_powers_0[i+1];
  }
  y_powers_0 = vector_powers(y, _num_multi_gate - 1);
  y_inv_powers_0 = vector_powers(y.inverse(), _num_multi_gate - 1);

  libff::enter_block("Verifier: flatten_wL, flatten_wR, flatten_wO");

  libff::enter_block("Verifier: flatten_wL, flatten_wR, flatten_wO only");
  flatten_constraints(_constraints, z, _num_multi_gate, flatten_wc, flatten_wL, flatten_wR, flatten_wO);

  libff::leave_block("Verifier: flatten_wL, flatten_wR, flatten_wO only");

  vector<Fr> l_x_3;
  vector<Fr> r_x_3;
  l_x_3 = hadmard_product(y_inv_powers_0, flatten_wR);
  r_x_3 = flatten_wL;
  delta = inner_product(l_x_3, r_x_3);

  // select random x
  transcript_g1.push_back(pi._pi_tilde);
  transcript_gt.push_back(pi._T_1);
  transcript_gt.push_back(pi._T_3);
  transcript_gt.push_back(pi._T_4);
  transcript_gt.push_back(pi._T_5);
  transcript_gt.push_back(pi._T_6);
  x = gen_random_field_elements_from_transcripts(transcript_g1, transcript_g2, transcript_gt, 1)[0];

  // select random x_t_tilde
  transcript_fr.push_back(pi._theta_x);
  transcript_fr.push_back(pi._mu);
  transcript_fr.push_back(pi._t_tilde);
  x_t_tilde = gen_random_field_elements_from_transcripts(transcript_fr, transcript_g1, transcript_g2, transcript_gt, (size_t)1)[0];

  // another pre computations
  vector<Fr> bold_wL, bold_wR, bold_wO;
  G1 WR_apos, WLO_apos;
  G1 P_G1, P_G1_apos;
  G1 P_G1_v, P_G1_apos_v;

  bold_wL = hadmard_product(y_inv_powers_0, flatten_wL);
  bold_wR = hadmard_product(y_inv_powers_0, flatten_wR);
  bold_wO = hadmard_product(y_inv_powers_0, flatten_wO);
  libff::leave_block("Verifier: flatten_wL, flatten_wR, flatten_wO");

  libff::enter_block("Verifier: wL, wR, wO");
  WR_apos = G1::zero();
  WLO_apos = G1::zero();
  for(size_t i = 1; i <= _num_multi_gate; i++) {
    WR_apos = WR_apos + bold_wR[i-1] * _srs.get_g1_beta_exp(i);
    WLO_apos = WLO_apos + (bold_wL[i-1] * x + bold_wO[i-1]) * _srs.get_g1_beta_exp((_N+1)+i);
  }

  std::vector<G1> g1tuple;
  std::vector<G2> g2tuple;
  CKzgKey kzg_ck;

  g1tuple = {_srs.get_g1_beta_exp(0)};
  g2tuple = {_srs.get_g2_beta_exp(0), _srs.get_g2_beta_exp(1)};
  kzg_ck = CKzgKey(g1tuple, g2tuple);

  assert(WR_apos == pi._commit_fwr_apos);
  result = kzg_vfyeval(kzg_ck, pi._commit_fwr_apos, pi._pi_kzg_fwr_apos);

  assert(WLO_apos == pi._commit_fwlo_apos);
  result &= kzg_vfyeval(kzg_ck, pi._commit_fwlo_apos, pi._pi_kzg_fwlo_apos);

  P_G1_v = G1::zero();
  P_G1_v = P_G1_v + (x ^ 1) * (pi._AI);
  P_G1_v = P_G1_v + (x ^ 1) * WR_apos;
  P_G1_v = P_G1_v + (x ^ 2) * pi._AO;
  P_G1_v = P_G1_v + (x ^ 3) * pi._K;
  P_G1_v = P_G1_v + WLO_apos;
  P_G1_v = P_G1_v + (-pi._mu) * _srs.get_g1_beta_exp(_N);
  P_G1_v = P_G1_v + _beta_powers_sum_neg;

  libff::leave_block("Verifier: wL, wR, wO");
  
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

  check_r_1 = _srs.gt ^ (pi._t_tilde - (x^2) * (delta + flatten_wc));
  check_r_2 = beta_N_plus_2_GT ^ (pi._theta_x);
  check_r_3 = ReducedPairing((x^2) * pi._pi_tilde, _srs.get_g2_beta_exp(0));

  check_l = check_l_1 * check_l_2;
  check_r = check_r_1 * check_r_2 * check_r_3;

  result &= (check_l == check_r);

#ifdef CRED_DEBUG
  libff::leave_block("/** varient of pointproof verify */");
  assert(result);
  libff::enter_block("/** inner product check for range proof **/");
#endif

  agenda.create_item("/** inner product check for range proof **/");

  std::vector<G1> g_vec;          // [[β^1]_1, ..., [β^D]_1]
  std::vector<G1> h_vec;          // [[y^(1-1) * β^(N+1+1)]_1, ..., [[y^(1-i) * β^(N+1+i)]_1], ..., [y^(1-D) * β^(N+1+D)]_1]

  for(size_t i = 1; i <= this->_num_multi_gate; i++) {
    g_vec.push_back(this->_srs.get_g1_beta_exp(i));
    h_vec.push_back(y_inv_powers_0[i-1] * this->_srs.get_g1_beta_exp((this->_N+1)+i));
  }

  this->_ipa_prove_sys._u = x_t_tilde * _srs.get_g1_beta_exp(_N);
  this->_ipa_prove_sys._P = P_G1_v;
  this->_ipa_prove_sys._c = pi._t_tilde;
  // P_G1_apos = P_G1 + pi._t_tilde * this->_ipa_prove_sys._u;

  agenda.create_item("Verify_base::IpaVerify");
  result &= this->_ipa_prove_sys.IpaMultiExpVerify(pi._pi_ipa, g_vec, h_vec);

  G1 g, h;
  g = pi._commit_fg;
  h = pi._commit_fh;
  agenda.create_item("Verify_improved::IpaVerify");
  result &= this->_ipa_prove_sys.IpaMultiExpVerify(pi._pi_ipa, g, h);
  result &= kzg_vfyeval(kzg_ck, pi._commit_fg, pi._pi_kzg_fg);
  result &= kzg_vfyeval(kzg_ck, pi._commit_fh, pi._pi_kzg_fh);
  agenda.mark_item_end("Verify_improved::IpaVerify");
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