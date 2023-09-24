#include <algorithm>
#include <assert.h>
#include <map>
#include <iostream>
#include <iterator>
#include "libff/common/utils.hpp"
#include "libff/algebra/curves/public_params.hpp"
#include "structs.hpp"
#include "prover.hpp"

/*
  以向量为系数的多项式
  存储问题
  运算问题：内积, 哈德玛积, 数乘
  l(X) = (_a_L - z * 1^D) + _K_L * X
  r(X) = hadmade( y^D, (_a_R) + z * 1^D + _K_R * X) ) + sum_(j \in _I) z^(1+j) * d_j
       = hadmade( y^D, (_a_R) + z * 1^D) ) + sum_(j \in _I) z^(1+j) * d_j + 
          hadmade( y^D, _K_R * X ) * X
  t0 = <l(x) 常数项系数, r(x) 常数项系数>
  t1 = <l(X) 常数项系数, r(x) 一次项系数> + <l(X) 一次项系数, r(x) 常数项系数>
  t2 = <l(X) 一次项系数, r(x) 一次项系数>
*/

namespace cred {

/************************** Prover ***************************/
const std::vector<Fr>& Prover::get_dj(size_t j) {
  return _ds.at(j);
}

void Prover::commit_set() {
  _r = Fr::random_element();
  _C = _r * this->_srs.get_g1_beta_exp(_N);  // item with index n havs exponent n+1
  // _C = G1::zero();  // item with index n havs exponent n+1
  for(size_t i = 1; i <= _S.m; i++) {
    _C = _C + _S.get_set_value(i) * this->_srs.get_g1_beta_exp(i);
  }

  // std::vector<Fr> exps = vector_powers(_srs.beta_powers[1], _S.m, false);
  // assert(_C == inner_product(_S.set_values, exps) * G1::one());
}


Prover::Prover(const CredSRS srs, const SET S, const Ranges& ranges, const IPAProveSystem& ipa_sys) {
  _srs = CredSRS(srs);
  _S = SET(S);
  _ranges = Ranges(ranges);
  _ipa_prove_sys = ipa_sys;

  _N = _srs.N;
  _D = _ranges._total_bits;
  _part_sums = _ranges._bits_part_sums;
  
  cout << "_D =" << _D << endl; 

  std::ostream_iterator<size_t> outie(cout, ", ");
  cout << "part_sums = \n" << endl;
  for(auto item: _part_sums) {
    cout << "key: " << item.first << " , value: " << item.second << endl;
  }

  commit_set();
}

CredProof Prover::prove() {

  CredProof pi;
  cout << "enter prove" << endl;
  Fr two = Fr::one() + Fr::one();
  std::vector<Fr> zero_vec;
  std::vector<Fr> one_vec;
  std::vector<Fr> z_vec;
  std::vector<Fr> z_powers_0; // index 0 with z^0
  std::vector<Fr> y_powers_0; // index 0 with y^1
  std::ostream_iterator<Fr>   outFr(cout, ", ");
  std::ostream_iterator<bool> outBo(cout, ", ");

  /************************* ds ***********************/
  zero_vec = std::vector<Fr>(_D, Fr::zero());
  
  size_t base;
  std::vector<Fr> dj_bold;

  for(size_t j: _ranges.get_indexs()) {
    const Range& range = _ranges.get_range(j);

    dj_bold = zero_vec;

    base  = _part_sums[j];
    dj_bold[base+0] = Fr::one();
    for(size_t i = 1; i < range.get_bit_len(); i++) {
      dj_bold[base+i] = two * dj_bold[base+i-1];
    }

    _ds.insert(std::pair<size_t, std::vector<Fr>>(j, dj_bold));
  }
  cout << "end ds" << endl;

  /************************* aL, aR ***********************/
  std::vector<Fr> values;
  libff::bit_vector values_bits, values_bits_neg;

  // for(size_t i: _ranges.get_indexs()) {
  //   values.push_back(_S.get_set_value(i));
  // }
  // values_bits = libff::convert_field_element_vector_to_bit_vector(values);
  // _a_L = libff::convert_bit_vector_to_field_element_vector<Fr>(values_bits);
  for(size_t i: _ranges.get_indexs()) {
    libff::bit_vector value_bits = libff::convert_field_element_to_bit_vector(_S.get_set_value(i));
    values_bits.insert(values_bits.end(), value_bits.begin(), value_bits.begin()+_ranges.get_range(i).get_bit_len());
  }

  _a_L = libff::convert_bit_vector_to_field_element_vector<Fr>(values_bits);
  _a_R = vector_sub(_a_L, std::vector<Fr>(_a_L.size(), Fr::one()));
  cout << "end aL, aR" << endl;

  /************************* r1, r2, KL, KR ***********************/
  _r_1 = Fr::random_element();
  _r_2 = Fr::random_element();
  for(size_t i = 1; i <= _D; i++) {
    _K_L.emplace_back(Fr::random_element());
    _K_R.emplace_back(Fr::random_element());
  }
  cout << "end r1, r2, KL, KR" << endl;

  /************************* A, K ***********************/
  cout << "_N = " << _N << endl; 
  _A = _r_1 * this->_srs.get_g1_beta_exp(_N);
  _K = _r_2 * this->_srs.get_g1_beta_exp(_N);
  for(size_t i = 1; i <= _D; i++) {
    _A = _A + _a_L[i-1] * this->_srs.get_g1_beta_exp(i) + _a_R[i-1] * this->_srs.get_g1_beta_exp(_N+1+i);
    _K = _K + _K_L[i-1] * this->_srs.get_g1_beta_exp(i) + _K_R[i-1] * this->_srs.get_g1_beta_exp(_N+1+i);
  }
  pi._A = _A;
  pi._K = _K;

  cout << "end A, K" << endl;

  /************************* Round 2 ***********************/
  /************************* y,z ***********************/
  generate_random_y_z(pi._A, pi._K, _y, _z);
  z_powers_0 = vector_powers(_z, _S.m+2); // 1, z, z^2, ... , z^m
  y_powers_0 = vector_powers(_y, _D-1);

  /************************* l(X) ***********************/
  std::vector<Fr> l_zero_term_coef, l_one_term_coef;

  one_vec = std::vector<Fr>(_D, Fr::one());
  z_vec   = numerical_mult(_z, one_vec);
  l_one_term_coef  = std::vector<Fr>(_K_L);
  l_zero_term_coef = vector_sub(_a_L, z_vec);
  cout << "  end l(X)" << endl;

  /************************* r(X) ***********************/
  std::vector<Fr> sum_z_pows_dot_dj, r_zero_term_coef, r_one_term_coef;

  sum_z_pows_dot_dj = zero_vec;
  for(size_t j: _ranges.get_indexs()) {
    sum_z_pows_dot_dj = vector_add(sum_z_pows_dot_dj, numerical_mult(z_powers_0[j+1], get_dj(j)));
  }

  r_one_term_coef  = hadmard_product(y_powers_0, _K_R);
  r_zero_term_coef = vector_add(sum_z_pows_dot_dj, hadmard_product(y_powers_0, vector_add(_a_R, z_vec)));
  cout << "  end r(X)" << endl;

  /************************* t0, t1, t2 ***********************/
  _t_0  = inner_product(l_zero_term_coef, r_zero_term_coef);
  _t_1  = inner_product(l_zero_term_coef, r_one_term_coef);
  _t_1 += inner_product(l_one_term_coef, r_zero_term_coef);
  _t_2  = inner_product(l_one_term_coef, r_one_term_coef);

  /************************* pi_tilde ***********************/
  G1 point_proof;

  _pi_tilde = G1::zero();
  for(size_t j: _ranges.get_indexs()) {
    point_proof = this->_r * _srs.get_g1_beta_exp(2 * _N + 1 - j);
    for(size_t i = 1; i <= _S.m; i++) {
      if(i != j) {
        point_proof = point_proof + _S.get_set_value(i) * this->_srs.get_g1_beta_exp(_N + 1 - j + i);
      }
    }
    _pi_tilde = _pi_tilde + z_powers_0[j+1] * point_proof;
  }
  cout << "end pi_tilde" << endl;
  
  /************************* r_3, r_4 ***********************/
  _r_3 = Fr::random_element();
  _r_4 = Fr::random_element();

  /************************* T1, T2 ***********************/
  GT beta_N_plus_2_GT = ReducedPairing(this->_srs.get_g1_beta_exp(_N+2), G2::one());

  _T_1 = ((_srs.gt) ^ (_t_1)) * (beta_N_plus_2_GT ^ (_r_3));
  _T_2 = ((_srs.gt) ^ (_t_2)) * (beta_N_plus_2_GT ^ (_r_4));

  pi._pi_tilde = _pi_tilde;
  pi._T_1 = _T_1;
  pi._T_2 = _T_2;

  cout << "end of round 2" << endl;
  /************************* Round 3 ***********************/
  /************************* x ***********************/
  generate_random_x(pi._pi_tilde, pi._T_1, pi._T_2, _x);

  /************************* l_bold, r_bold, t_tilde ***********************/
  _t_tilde = _t_0 + _t_1 * _x + _t_2 * (_x^2);
  
  /************************* r_x, mu ***********************/
  _r_x  = _r_3 * _x + _r_4 * (_x^2);
  _mu   = _r_1      + _r_2 * _x;

  pi._t_tilde = _t_tilde;
  pi._r_x = _r_x;
  pi._mu = _mu;
  

  /************************* [P]_1 ***********************/
  std::vector<Fr> y_inv_powers_0;
  
  y_inv_powers_0 = vector_powers(_y.inverse(), _D - 1);
  _l_bold = vector_add(l_zero_term_coef, numerical_mult(_x, l_one_term_coef));
  _r_bold = vector_add(r_zero_term_coef, numerical_mult(_x, r_one_term_coef));

  std::vector<G1> g_vec, h_vec;
  for(size_t i = 0; i < this->_D; i++) {
    g_vec.push_back(this->_srs.get_g1_beta_exp(i+1));
    h_vec.push_back(y_inv_powers_0[i] * this->_srs.get_g1_beta_exp(this->_N+2+i));
  }

  // _P_G1 = G1::zero();
  _P_G1 = pi._t_tilde * this->_ipa_prove_sys._u;
  for(size_t i = 0; i < _D; i++) {
    _P_G1 = _P_G1 + _l_bold[i] * g_vec[i];
    _P_G1 = _P_G1 + _r_bold[i] * h_vec[i];
  }
  /************************* buller_ipa todo ***********************/
  this->_ipa_prove_sys._P = this->_P_G1;
  this->_ipa_prove_sys._c = pi._t_tilde;
  pi._pi_ipa = this->_ipa_prove_sys.IpaProve(g_vec, h_vec, _l_bold, _r_bold, pi._t_tilde);

  assert(t0_check(two, one_vec, z_powers_0, y_powers_0));
  assert(easy_check(z_powers_0));
  cout << "pass: assert(left == right_1 * right_2);" << endl;

  return pi;
}

bool Prover::t0_check(const Fr& two, const vector<Fr>& one_vec, const vector<Fr>& z_powers_0, const vector<Fr>& y_powers_0) {
  Fr sum_z_pows_dot_inner = Fr::zero();
  for(size_t j: _ranges.get_indexs()) {
    Fr two_pow            = two ^ (_ranges.get_range(j).get_bit_len());
    sum_z_pows_dot_inner += z_powers_0[j+2] * (two_pow - Fr::one());
  }

  // Fr t0_1, t0_2, t0_3;
  t0_1 = Fr::zero();
  for(size_t j: _ranges.get_indexs()) {
    t0_1 += z_powers_0[j+1] * inner_product(_a_L, get_dj(j)) ;
  }
  t0_2 = (z_powers_0[1]-z_powers_0[2]) * inner_product(one_vec, y_powers_0);
  t0_3 = sum_z_pows_dot_inner;

  return _t_0 == t0_1 + t0_2 - t0_3;
}

/*
  point proof 多元素聚合检查, 不考虑零知识性
*/
bool Prover::easy_check(const std::vector<Fr>& z_powers_0) {

  G2 sum_z_pows_beta_pows;
  sum_z_pows_beta_pows = G2::zero();
  for(size_t j: _ranges.get_indexs()) {
    sum_z_pows_beta_pows = sum_z_pows_beta_pows + z_powers_0[1+j] * this->_srs.get_g2_beta_exp(_N + 1 - j);
  }

  GT left = ReducedPairing(_C, sum_z_pows_beta_pows);
  GT right_1 = ReducedPairing(_pi_tilde, G2::one());
  Fr right_2_Fr = Fr::zero();
  for(size_t j: _ranges.get_indexs()) {
    right_2_Fr += _S.get_set_value(j) * z_powers_0[j+1];
  }
  GT right_2 = _srs.gt ^ right_2_Fr;

  return left == right_1 * right_2;
}

};