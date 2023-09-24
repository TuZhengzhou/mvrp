#include "structs.hpp"
#include "prover.hpp"
#include "verifier.hpp"

namespace cred {

/************************** Verifier ***************************/
Verifier::Verifier(const CredSRS srs, const G1 C, const size_t set_size, const Ranges& ranges, const IPAProveSystem& ipa_sys) {
  _srs = CredSRS(srs);
  _set_size = set_size;
  _ranges = Ranges(ranges);
  _ipa_prove_sys = ipa_sys;

  _N = _srs.N;
  _D = _ranges._total_bits;
  _part_sums = _ranges._bits_part_sums;

  _C = G1(C);
}

bool Verifier::verify(const CredProof& pi, const Prover& prover) {

  libff::enter_block("Verifier::verify");

  Fr two = Fr::one() + Fr::one();
  Fr two_pow;
  std::vector<Fr> z_powers_0; // index 0 with z^0
  std::vector<Fr> one_vec;
  std::vector<Fr> y_powers_0; // index 0 with y^0

  /************************* y,z ***********************/
  generate_random_y_z(pi._A, pi._K, _y, _z);
  z_powers_0 = vector_powers(_z, _set_size+2); // 1, z, z^2, ... , z^m

  /************************* x ***********************/
  generate_random_x(pi._pi_tilde, pi._T_1, pi._T_2, _x);

  /************************* delta ***********************/
  one_vec = std::vector<Fr>(_D, Fr::one());
  y_powers_0 = vector_powers(_y, _D-1);

  Fr sum_z_pows_dot_inner = Fr::zero();
  for(size_t j: _ranges.get_indexs()) {
    two_pow = two ^ (_ranges.get_range(j).get_bit_len());
    sum_z_pows_dot_inner       += z_powers_0[j+2] * (two_pow - Fr::one());
  }
  _delta = (z_powers_0[1] - z_powers_0[2]) * inner_product(one_vec, y_powers_0) - sum_z_pows_dot_inner;

  /************************* left items ***********************/
  // GT _left_1, _left_2, _left_3, _left;
  _left_1 = pi._T_1 ^ _x;
  _left_2 = pi._T_2 ^ (_x ^ 2);

  G2 sum_z_pows_beta_pows;
  sum_z_pows_beta_pows = G2::zero();
  for(size_t j: _ranges.get_indexs()) {
    sum_z_pows_beta_pows = sum_z_pows_beta_pows + z_powers_0[1+j] * this->_srs.get_g2_beta_exp(_N + 1 - j);
  }
  _left_3 = ReducedPairing(_C, sum_z_pows_beta_pows);
  _left = _left_1 * _left_2 * _left_3;

  /************************* right items ***********************/
  // GT _right_1, _right_2, _right_3, _right;
  GT beta_N_plus_2_GT = ReducedPairing(this->_srs.get_g1_beta_exp(_N+2), G2::one());

  _right_1 = _srs.gt ^ (pi._t_tilde - _delta);
  _right_2 = beta_N_plus_2_GT ^ (pi._r_x);
  _right_3 = ReducedPairing(pi._pi_tilde, G2::one());
  _right = _right_1 * _right_2 * _right_3;

  assert(_left == _right);

  /************************* [F] ***********************/
  size_t bit_len, base;
  std::vector<Fr> y_inv_powers_0;
  std::vector<G1> g_vec, h_vec;
  std::vector<Fr> l_apos, r_apos; // - z \cdot 1^D, y_bold ^ D \circuit (z\cdot 1^d) + \sum ...
  std::vector<Fr> sum_z_pow_dj;
  G1 F, F1, P_apos, inner_sum;

  y_inv_powers_0 = vector_powers(_y.inverse(), _D - 1);

  for(size_t i = 0; i < this->_D; i++) {
    g_vec.push_back(this->_srs.get_g1_beta_exp(1+i));
    h_vec.push_back(y_inv_powers_0[i] * this->_srs.get_g1_beta_exp(this->_N+2+i));
  }

  l_apos = std::vector<Fr>(_D, - _z * Fr::one());
  r_apos = hadmard_product(y_powers_0, std::vector<Fr>(_D, _z * Fr::one()));
  sum_z_pow_dj = std::vector<Fr>(_D, Fr::one());

  base = 0;
  for(auto j: this->_ranges.get_indexs()) {

    bit_len = this->_ranges.get_range(j).get_bit_len();
    sum_z_pow_dj[base] = z_powers_0[j + 1];

    for(size_t i = 1; i < bit_len; i++) {
      sum_z_pow_dj[base + i] = sum_z_pow_dj[base + i - 1] * two;
    }
    base += bit_len;
  }

  F = G1::zero();
  for(size_t i = 0; i < _D; i++) {
    F = F + l_apos[i] * g_vec[i] + (r_apos[i] + sum_z_pow_dj[i]) * h_vec[i];
  }

  P_apos = pi._A + (_x * pi._K) + F + (-pi._mu) * _srs.get_g1_beta_exp(_N);

  this->_ipa_prove_sys._P = P_apos;
  this->_ipa_prove_sys._c = pi._t_tilde;
  assert(this->_ipa_prove_sys.IpaVerify(pi._pi_ipa, g_vec, h_vec));

  libff::leave_block("Verifier::verify");
  return _left == _right;
}


};
