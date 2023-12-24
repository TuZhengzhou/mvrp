#include <algorithm>
#include <assert.h>
#include <map>
#include <iostream>
#include <iterator>
#include "structs.hpp"
#include "range_prover.hpp"


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

void Prover::commit_set() {
  _r = Fr::random_element();
  _C = _r * this->_srs.get_g1_beta_exp(_N);  // item with index n havs exponent n+1
  for(size_t i = 1; i <= _S.m; i++) {
    _C = _C + _S.get_set_value(i) * this->_srs.get_g1_beta_exp(i);
  }
}

void Prover::point_proof_pre_compute() {

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

const G1& Prover::point_proof(size_t idx) {
#ifdef CRED_DEBUG
  assert(idx >= 1 && idx <= _S.set_values.size());
#endif
  return _point_proofs[idx-1];
}



Prover::Prover(const CCredSRS& srs, const CSet& S, const CRanges& ranges, const IPAProveSystem& ipa_sys) {
  _srs = CCredSRS(srs);
  _S = CSet(S);
  _ranges = CRanges(ranges);
  _ipa_prove_sys = ipa_sys;

  _N = _srs.N;
  _D = _ranges._total_bits;
  _part_sums = _ranges._bits_part_sums;

#ifdef CRED_DEBUG
  cout << "_D =" << _D << endl; 
  std::ostream_iterator<size_t> outie(cout, ", ");
  for(auto item: _part_sums) {
    cout << "part_sums mapping [key: " << item.first << " , value: " << item.second << "]" << endl;
  }
#endif

  commit_set();
  point_proof_pre_compute();
}

CRangeProof Prover::prove(Agenda& agenda, const bool improved) {
  if(improved) 
    return prove_improved(agenda);
  return prove_base(agenda);
}

CRangeProof Prover::prove_base(Agenda& agenda) {

  agenda.create_item("Prove_base");
  
  CRangeProof pi; // global proof pi
  std::ostream_iterator<Fr>   outFr(cout, ", ");
  std::ostream_iterator<bool> outBo(cout, ", ");

  vector<Fr> transcript_fr;
  vector<G1> transcript_g1;
  vector<G2> transcript_g2;
  vector<GT> transcript_gt;
  vector<Fr> random_challenges;

  // round 1
  Fr _r_1, _r_2;
  std::vector<Fr> _a_L, _a_R, _K_L, _K_R;
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
  
  agenda.create_item("bullet::Round1");
  agenda.create_item("Prove_base::Round1");
  /**
    Round1: commit_for_ipa, output [A]_1, [K]_1
  **/
  size_t base, i;
  Fr two;                   // 2 in Fr
  std::vector<Fr> zero_vec; // [0, ..., 0]
  std::vector<Fr> dj_bold;  // [1, 2, ..., 2^bitlen_i] for some value_i in set, bitlen_i is the bitlen of value_i  

  two = Fr::one() + Fr::one();
  zero_vec = std::vector<Fr>(_D, Fr::zero());

  const bool round_1 = true;
  if(round_1) {

    /** for each j in I, compute vector dj (Prover._ds) **/
#ifdef CRED_DEBUG
    libff::enter_block("/** for each j in I, compute vector dj (Prover._ds) **/");
#endif

#ifdef CRED_DEBUG
    libff::leave_block("/** for each j in I, compute vector dj (Prover._ds) **/");
#endif



    /** compute a_L \in {0,1}^D and a_R \in Z^D_p satisfies <a_L,dj> = vj and a_R = a_L - 1^D **/
#ifdef CRED_DEBUG
    libff::enter_block("/** compute a_L in {0,1}^D and a_R in Z^D_p satisfies <a_L,dj> = vj and a_R = a_L - 1^D **/");
#endif

    libff::bit_vector value_bits;       // value_i into bits
    libff::bit_vector values_bits;      // [value_1, ..., value_n] into bits(total D bits)
    libff::bit_vector values_bits_neg;  // values_bits - 1^D

    for(auto i: _ranges.get_indexs()) {
      value_bits = libff::convert_field_element_to_bit_vector(_S.get_set_value(i));
      values_bits.insert(values_bits.end(), value_bits.begin(), value_bits.begin()+_ranges.get_range(i).get_bit_len());
    }

    _a_L = libff::convert_bit_vector_to_field_element_vector<Fr>(values_bits);
    _a_R = vector_sub(_a_L, std::vector<Fr>(_a_L.size(), Fr::one()));

#ifdef CRED_DEBUG
    libff::leave_block("/** compute a_L in {0,1}^D and a_R in Z^D_p satisfies <a_L,dj> = vj and a_R = a_L - 1^D **/");
#endif

    /** sample r1, r2 ←$ Zp, KL, KR ←$ Z^D_p **/
    _r_1 = Fr::random_element();
    _r_2 = Fr::random_element();
    _K_L = std::vector<Fr>(_D);
    _K_R = std::vector<Fr>(_D);
    // Fr::random_e
    for(i = 0; i < _D; i++) {
      _K_L[i] = Fr::random_element();
      _K_R[i] = Fr::random_element();
    }

    /** compute 
     * [A]1 = [D∑i=1(aL[i]β^i + aR[i]β^(N+1+i)) + r_1 β^N ]_1 and 
     * [K]1 = [D∑i=1(KL[i]β^i + KR[i]β^(N+1+i)) + r_2 β^N ]_1
    **/
#ifdef CRED_DEBUG
    libff::enter_block("/** compute [A]_1, [K]_1 **/");
#endif

    _A = _r_1 * this->_srs.get_g1_beta_exp(_N);
    _K = _r_2 * this->_srs.get_g1_beta_exp(_N);
    for(size_t i = 1; i <= _D; i++) {
      _A = _A + _a_L[i-1] * this->_srs.get_g1_beta_exp(i) + _a_R[i-1] * this->_srs.get_g1_beta_exp(_N+1+i);
      _K = _K + _K_L[i-1] * this->_srs.get_g1_beta_exp(i) + _K_R[i-1] * this->_srs.get_g1_beta_exp(_N+1+i);
    }

#ifdef CRED_DEBUG
    libff::leave_block("/** compute [A]_1, [K]_1 **/");
#endif

    pi._A = _A;
    pi._K = _K;
  }
  agenda.mark_item_end("bullet::Round1");
  agenda.mark_item_end("Prove_base::Round1");



  agenda.create_item("bullet::Round2_plus_pi_tilde");
  agenda.create_item("Prove_base::Round2");
  /**
    Round2: point_prove, output [pi_tilde]_1, [T_1]_T, [T_2]_T
  **/
  std::vector<Fr> one_vec;    // [1, ..., 1]
  std::vector<Fr> z_vec;      // [z, ..., z]
  std::vector<Fr> z_powers_0; // index 0 with z^0, [1, z, z^2, ... , z^m, z^(m+1)]
  std::vector<Fr> y_powers_0; // index 0 with y^1, [1, y, y^2, ..., y^(D-1)]

  std::vector<Fr> sum_z_pows_cdot_dj;                 // ∑ j in I {z^(1+j) \cdot d_j}
  std::vector<Fr> l_zero_term_coef, l_one_term_coef;  // coefs of l(x)
  std::vector<Fr> r_zero_term_coef, r_one_term_coef;  // coefs of r(x)
  GT beta_N_plus_2_GT;

  // G1 point_proof;   // point proof of one item
  
  const bool round_2 = true;
  if(round_2) {

#ifdef CRED_DEBUG
    libff::enter_block("/** Round2: point_prove, output [pi_tilde]_1, [T_1]_T, [T_2]_T **/");
#endif

    /**  y, z ←$ hash([A]_1, [K]_1) **/
    transcript_g1.push_back(pi._A);
    transcript_g1.push_back(pi._K);
    random_challenges = gen_random_field_elements_from_transcripts(transcript_g1, transcript_g2, transcript_gt, 2);
    _y = random_challenges[0];
    _z = random_challenges[1];

    z_powers_0 = vector_powers(_z, _S.m+2); // 1, z, z^2, ... , z^m, z^(m+1)
    y_powers_0 = vector_powers(_y, _D-1);

    one_vec = std::vector<Fr>(_D, Fr::one());
    z_vec   = numerical_mult(_z, one_vec);

    /** compute coefs of l(x) **/
    l_one_term_coef  = std::vector<Fr>(_K_L);
    l_zero_term_coef = vector_sub(_a_L, z_vec);

    /** compute coefs of r(x) **/
    libff::enter_block("/** compute sum_z_pows_cdot_dj **/");
    sum_z_pows_cdot_dj = zero_vec;
    for(auto j: _ranges.get_indexs()) {
      // sum_z_pows_cdot_dj = vector_add(sum_z_pows_cdot_dj, numerical_mult(z_powers_0[j+1], get_dj(j)));
      base = _ranges.get_part_sum(j);
      auto bit_len = _ranges.get_range(j).get_bit_len();
      auto two_exp = Fr::one();
      for(i = 0; i < bit_len; i++) {
        sum_z_pows_cdot_dj[base+i] += z_powers_0[j+1] * two_exp;
        two_exp *= two;
      }
    }
    libff::leave_block("/** compute sum_z_pows_cdot_dj **/");
    r_one_term_coef  = hadmard_product(y_powers_0, _K_R);
    r_zero_term_coef = vector_add(sum_z_pows_cdot_dj, hadmard_product(y_powers_0, vector_add(_a_R, z_vec)));

    /** t(x) = l(x) * r(x), compute coefs of t(x) **/
    _t_0  = inner_product(l_zero_term_coef, r_zero_term_coef);
    _t_1  = inner_product(l_zero_term_coef, r_one_term_coef );
    _t_1 += inner_product(l_one_term_coef , r_zero_term_coef);
    _t_2  = inner_product(l_one_term_coef , r_one_term_coef );

    /** compute pi_tilde: pointproof of multi items **/
    agenda.create_item("Prove::_pi_tilde");
    _pi_tilde = G1::zero();
    for(auto j: _ranges.get_indexs()) {
      _pi_tilde = _pi_tilde + z_powers_0[j+1] * this->point_proof(j);
    }
    agenda.mark_item_end("Prove::_pi_tilde");
    
    /************************* r_3, r_4 ***********************/
    _r_3 = Fr::random_element();
    _r_4 = Fr::random_element();

    /************************* T1, T2 ***********************/
    // beta_N_plus_2_GT = ReducedPairing(this->_srs.get_g1_beta_exp(_N+2), _srs.get_g2_beta_exp(0)));
  beta_N_plus_2_GT = ReducedPairing(this->_srs.get_g1_beta_exp(_N), this->_srs.get_g2_beta_exp(2));

    _T_1 = ((_srs.gt) ^ (_t_1)) * (beta_N_plus_2_GT ^ (_r_3));
    _T_2 = ((_srs.gt) ^ (_t_2)) * (beta_N_plus_2_GT ^ (_r_4));

#ifdef CRED_DEBUG
    libff::leave_block("/** Round2: point_prove, output [pi_tilde]_1, [T_1]_T, [T_2]_T **/");
#endif

    pi._pi_tilde = _pi_tilde;
    pi._T_1 = _T_1;
    pi._T_2 = _T_2;
  }
  agenda.mark_item_end("bullet::Round2_plus_pi_tilde");
  agenda.mark_item_end("Prove_base::Round2");
  

  agenda.create_item("bullet::Round3");
  agenda.create_item("Prove_base::Round3");
  /**
    Round3: range_prove, output \mu, r_x, pi_ipa, t_tilde = t(x)
  **/
  std::vector<Fr> y_inv_powers_0; // [1, y^-1, ..., y^-(D-1)]
  std::vector<G1> g_vec;          // [[β^1]_1, ..., [β^D]_1]
  std::vector<G1> h_vec;          // [[y^(1-1) * β^(N+1+1)]_1, ..., [[y^(1-i) * β^(N+1+i)]_1], ..., [y^(1-D) * β^(N+1+D)]_1]

  y_inv_powers_0 = vector_powers(_y.inverse(), _D - 1);
  for(size_t i = 0; i < this->_D; i++) {
    g_vec.push_back(this->_srs.get_g1_beta_exp(i+1));
    h_vec.push_back(y_inv_powers_0[i] * this->_srs.get_g1_beta_exp(this->_N+2+i));
  }

  const bool round_3 = true;
  if(round_3) {

#ifdef CRED_DEBUG
    libff::enter_block("/** Round3: range_prove, output \\mu, r_x, pi_ipa, t_tilde = t(x) **/");
#endif

    /** x ←$ hash([pi_tilde]_1, [T1]_T, [T2]_T) **/
    transcript_g1.push_back(pi._pi_tilde);
    transcript_gt.push_back(pi._T_1);
    transcript_gt.push_back(pi._T_2);
    _x = gen_random_field_elements_from_transcripts(transcript_g1, transcript_g2, transcript_gt, 1)[0];

    /** compute l(x), r(x), and t(x) **/
    _l_bold = vector_add(l_zero_term_coef, numerical_mult(_x, l_one_term_coef));
    _r_bold = vector_add(r_zero_term_coef, numerical_mult(_x, r_one_term_coef));
    _t_tilde = _t_0 + _t_1 * _x + _t_2 * (_x^2);
#ifdef CRED_DEBUG
    assert(_t_tilde == inner_product(_l_bold, _r_bold));
#endif

    /** compute \mu = r1 + r2 * x, rx = r3 * x + r4 * x^2 **/
    _mu   = _r_1      + _r_2 * _x;
    _r_x  = _r_3 * _x + _r_4 * (_x^2);

    /** compute  [P]1 = [D∑i=1( l[i] * β^i + r[i] * y^(−i+1) * β^(N+1+i) )] **/

    _P_G1 = _t_tilde * this->_ipa_prove_sys._u; // ???
    for(size_t i = 0; i < _D; i++) {
      _P_G1 = _P_G1 + _l_bold[i] * g_vec[i];
      _P_G1 = _P_G1 + _r_bold[i] * h_vec[i];
    }
    /************************* buller_ipa ***********************/
    this->_ipa_prove_sys._P = _P_G1;
    this->_ipa_prove_sys._c = _t_tilde;
    agenda.create_item("Prove::IpaProve");
    auto pi_ipa             = this->_ipa_prove_sys.IpaProve(g_vec, h_vec, _l_bold, _r_bold, _t_tilde);
    agenda.mark_item_end("Prove::IpaProve");

#ifdef CRED_DEBUG
    libff::leave_block("/** Round3: range_prove, output \\mu, r_x, pi_ipa, t_tilde = t(x) **/");
#endif

    pi._mu = _mu;
    pi._r_x = _r_x;
    pi._t_tilde = _t_tilde;
    pi._pi_ipa = pi_ipa;
  }

  agenda.mark_item_end("bullet::Round3");
  agenda.mark_item_end("Prove_base::Round3");

  agenda.mark_item_end("Prove_base");
  agenda.mark_proof_size(pi);
  agenda.mark_ipa_size(pi);

  agenda.mark_mem_usage("Prove_base");

  return pi;
};

CRangeProof Prover::prove_improved(Agenda& agenda) {

  libff::enter_block("Prover::prove_improved()");
  agenda.create_item("Prove_improved");

  CRangeProof pi; // global proof pi
  std::ostream_iterator<Fr>   outFr(cout, ", ");
  std::ostream_iterator<bool> outBo(cout, ", ");

  vector<Fr> transcript_fr;
  vector<G1> transcript_g1;
  vector<G2> transcript_g2;
  vector<GT> transcript_gt;
  vector<Fr> random_challenges;

  // round 1
  Fr _r_1, _r_2;
  std::vector<Fr> _a_L, _a_R, _K_L, _K_R;
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

  agenda.create_item("Prove_improved::Round1");
  /**
    Round1: commit_for_ipa, output [A]_1, [K]_1
  **/
  size_t base, i;
  Fr two;                   // 2 in Fr
  std::vector<Fr> zero_vec; // [0, ..., 0]
  std::vector<Fr> dj_bold;  // [1, 2, ..., 2^bitlen_i] for some value_i in set, bitlen_i is the bitlen of value_i  

  two = Fr::one() + Fr::one();
  zero_vec = std::vector<Fr>(_D, Fr::zero());

  const bool round_1 = true;
  if(round_1) {

    /** for each j in I, compute vector dj (Prover._ds) **/
#ifdef CRED_DEBUG
    libff::enter_block("/** for each j in I, compute vector dj (Prover._ds) **/");
#endif

#ifdef CRED_DEBUG
    libff::leave_block("/** for each j in I, compute vector dj (Prover._ds) **/");
#endif



    /** compute a_L \in {0,1}^D and a_R \in Z^D_p satisfies <a_L,dj> = vj and a_R = a_L - 1^D **/
#ifdef CRED_DEBUG
    libff::enter_block("/** compute a_L in {0,1}^D and a_R in Z^D_p satisfies <a_L,dj> = vj and a_R = a_L - 1^D **/");
#endif

    libff::bit_vector value_bits;       // value_i into bits
    libff::bit_vector values_bits;      // [value_1, ..., value_n] into bits(total D bits)
    libff::bit_vector values_bits_neg;  // values_bits - 1^D

    for(auto i: _ranges.get_indexs()) {
      value_bits = libff::convert_field_element_to_bit_vector(_S.get_set_value(i));
      values_bits.insert(values_bits.end(), value_bits.begin(), value_bits.begin()+_ranges.get_range(i).get_bit_len());
    }

    _a_L = libff::convert_bit_vector_to_field_element_vector<Fr>(values_bits);
    _a_R = vector_sub(_a_L, std::vector<Fr>(_a_L.size(), Fr::one()));

#ifdef CRED_DEBUG
    libff::leave_block("/** compute a_L in {0,1}^D and a_R in Z^D_p satisfies <a_L,dj> = vj and a_R = a_L - 1^D **/");
#endif

    /** sample r1, r2 ←$ Zp, KL, KR ←$ Z^D_p **/
    _r_1 = Fr::random_element();
    _r_2 = Fr::random_element();
    _K_L = std::vector<Fr>(_D);
    _K_R = std::vector<Fr>(_D);
    for(i = 0; i < _D; i++) {
      _K_L[i] = Fr::random_element();
      _K_R[i] = Fr::random_element();
    }

    /** compute 
     * [A]1 = [D∑i=1(aL[i]β^i + aR[i]β^(N+1+i)) + r_1 β^N ]_1 and 
     * [K]1 = [D∑i=1(KL[i]β^i + KR[i]β^(N+1+i)) + r_2 β^N ]_1
    **/
#ifdef CRED_DEBUG
    libff::enter_block("/** compute [A]_1, [K]_1 **/");
#endif

    _A = _r_1 * this->_srs.get_g1_beta_exp(_N);
    _K = _r_2 * this->_srs.get_g1_beta_exp(_N);
    for(size_t i = 1; i <= _D; i++) {
      _A = _A + _a_L[i-1] * this->_srs.get_g1_beta_exp(i) + _a_R[i-1] * this->_srs.get_g1_beta_exp(_N+1+i);
      _K = _K + _K_L[i-1] * this->_srs.get_g1_beta_exp(i) + _K_R[i-1] * this->_srs.get_g1_beta_exp(_N+1+i);
    }

#ifdef CRED_DEBUG
    libff::leave_block("/** compute [A]_1, [K]_1 **/");
#endif

    pi._A = _A;
    pi._K = _K;
  }
  agenda.mark_item_end("Prove_improved::Round1");



  agenda.create_item("Prove_improved::Round2");
  /**
    Round2: point_prove, output [pi_tilde]_1, [T_1]_T, [T_2]_T
  **/
  std::vector<Fr> one_vec;    // [1, ..., 1]
  std::vector<Fr> z_vec;      // [z, ..., z]
  std::vector<Fr> z_powers_0; // index 0 with z^0, [1, z, z^2, ... , z^m, z^(m+1)]
  std::vector<Fr> y_powers_0; // index 0 with y^1, [1, y, y^2, ..., y^(D-1)]

  std::vector<Fr> sum_z_pows_cdot_dj;                 // ∑ j in I {z^(1+j) \cdot d_j}
  std::vector<Fr> l_zero_term_coef, l_one_term_coef;  // coefs of l(x)
  std::vector<Fr> r_zero_term_coef, r_one_term_coef;  // coefs of r(x)
  GT beta_N_plus_2_GT;
  
  const bool round_2 = true;
  if(round_2) {

#ifdef CRED_DEBUG
    libff::enter_block("/** Round2: point_prove, output [pi_tilde]_1, [T_1]_T, [T_2]_T **/");
#endif

    /**  y, z ←$ hash([A]_1, [K]_1) **/
    transcript_g1.push_back(pi._A);
    transcript_g1.push_back(pi._K);
    random_challenges = gen_random_field_elements_from_transcripts(transcript_g1, transcript_g2, transcript_gt, 2);
    _y = random_challenges[0];
    _z = random_challenges[1];

    z_powers_0 = vector_powers(_z, _S.m+2); // 1, z, z^2, ... , z^m, z^(m+1)
    y_powers_0 = vector_powers(_y, _D-1);

    one_vec = std::vector<Fr>(_D, Fr::one());
    z_vec   = numerical_mult(_z, one_vec);

    /** compute coefs of l(x) **/
    l_one_term_coef  = std::vector<Fr>(_K_L);
    l_zero_term_coef = vector_sub(_a_L, z_vec);

    /** compute coefs of r(x) **/
    libff::enter_block("/** compute sum_z_pows_cdot_dj **/");
    sum_z_pows_cdot_dj = zero_vec;
    for(auto j: _ranges.get_indexs()) {
      base = _ranges.get_part_sum(j);
      auto bit_len = _ranges.get_range(j).get_bit_len();
      auto two_exp = Fr::one();
      for(i = 0; i < bit_len; i++) {
        sum_z_pows_cdot_dj[base+i] += z_powers_0[j+1] * two_exp;
        two_exp *= two;
      }
    }
    libff::leave_block("/** compute sum_z_pows_cdot_dj **/");
    r_one_term_coef  = hadmard_product(y_powers_0, _K_R);
    r_zero_term_coef = vector_add(sum_z_pows_cdot_dj, hadmard_product(y_powers_0, vector_add(_a_R, z_vec)));

    /** t(x) = l(x) * r(x), compute coefs of t(x) **/
    _t_0  = inner_product(l_zero_term_coef, r_zero_term_coef);
    _t_1  = inner_product(l_zero_term_coef, r_one_term_coef );
    _t_1 += inner_product(l_one_term_coef , r_zero_term_coef);
    _t_2  = inner_product(l_one_term_coef , r_one_term_coef );

    /** compute pi_tilde: pointproof of multi items **/
    _pi_tilde = G1::zero();
    for(auto j: _ranges.get_indexs()) {
      _pi_tilde = _pi_tilde + z_powers_0[j+1] * this->point_proof(j);
    }
    
    /************************* r_3, r_4 ***********************/
    _r_3 = Fr::random_element();
    _r_4 = Fr::random_element();

    /************************* T1, T2 ***********************/
  beta_N_plus_2_GT = ReducedPairing(this->_srs.get_g1_beta_exp(_N), this->_srs.get_g2_beta_exp(2));

    _T_1 = ((_srs.gt) ^ (_t_1)) * (beta_N_plus_2_GT ^ (_r_3));
    _T_2 = ((_srs.gt) ^ (_t_2)) * (beta_N_plus_2_GT ^ (_r_4));

#ifdef CRED_DEBUG
    libff::leave_block("/** Round2: point_prove, output [pi_tilde]_1, [T_1]_T, [T_2]_T **/");
#endif

    pi._pi_tilde = _pi_tilde;
    pi._T_1 = _T_1;
    pi._T_2 = _T_2;
  }
  agenda.mark_item_end("Prove_improved::Round2");



  agenda.create_item("Prove_improved::Round3");
  /**
    Round3: range_prove, output \mu, r_x, pi_ipa, t_tilde = t(x)
  **/
  std::vector<Fr> y_inv_powers_0; // [1, y^-1, ..., y^-(D-1)]
  std::vector<G1> g_vec;          // [[β^1]_1, ..., [β^D]_1]
  std::vector<G1> h_vec;          // [[y^(1-1) * β^(N+1+1)]_1, ..., [[y^(1-i) * β^(N+1+i)]_1], ..., [y^(1-D) * β^(N+1+D)]_1]

  y_inv_powers_0 = vector_powers(_y.inverse(), _D - 1);
  for(size_t i = 0; i < this->_D; i++) {
    g_vec.push_back(this->_srs.get_g1_beta_exp(i+1));
    h_vec.push_back(y_inv_powers_0[i] * this->_srs.get_g1_beta_exp(this->_N+2+i));
  }

  const bool round_3 = true;
  if(round_3) {

#ifdef CRED_DEBUG
    libff::enter_block("/** Round3: range_prove, output \\mu, r_x, pi_ipa, t_tilde = t(x) **/");
#endif

    /** x ←$ hash([pi_tilde]_1, [T1]_T, [T2]_T) **/
    transcript_g1.push_back(pi._pi_tilde);
    transcript_gt.push_back(pi._T_1);
    transcript_gt.push_back(pi._T_2);
    _x = gen_random_field_elements_from_transcripts(transcript_g1, transcript_g2, transcript_gt, 1)[0];

    /** compute l(x), r(x), and t(x) **/
    _l_bold = vector_add(l_zero_term_coef, numerical_mult(_x, l_one_term_coef));
    _r_bold = vector_add(r_zero_term_coef, numerical_mult(_x, r_one_term_coef));
    _t_tilde = _t_0 + _t_1 * _x + _t_2 * (_x^2);
#ifdef CRED_DEBUG
    assert(_t_tilde == inner_product(_l_bold, _r_bold));
#endif

    /** compute \mu = r1 + r2 * x, rx = r3 * x + r4 * x^2 **/
    _mu   = _r_1      + _r_2 * _x;
    _r_x  = _r_3 * _x + _r_4 * (_x^2);

    /** compute  [P]1 = [D∑i=1( l[i] * β^i + r[i] * y^(−i+1) * β^(N+1+i) )] **/

    _P_G1 = _t_tilde * this->_ipa_prove_sys._u; // ???
    for(size_t i = 0; i < _D; i++) {
      _P_G1 = _P_G1 + _l_bold[i] * g_vec[i];
      _P_G1 = _P_G1 + _r_bold[i] * h_vec[i];
    }
    /************************* buller_ipa ***********************/
    this->_ipa_prove_sys._P = _P_G1;
    this->_ipa_prove_sys._c = _t_tilde;

    agenda.create_item("Prove::IpaProve");
    auto pi_ipa             = this->_ipa_prove_sys.IpaProve(g_vec, h_vec, _l_bold, _r_bold, _t_tilde);
    agenda.mark_item_end("Prove::IpaProve");

#ifdef CRED_DEBUG
    libff::leave_block("/** Round3: range_prove, output \\mu, r_x, pi_ipa, t_tilde = t(x) **/");
#endif

    pi._mu = _mu;
    pi._r_x = _r_x;
    pi._t_tilde = _t_tilde;
    pi._pi_ipa = pi_ipa;
  }






  /** 
   * kzg proof for polynomial fzy 
  */
  /** compute polynomial fzy */
  size_t bit_len;
  std::vector<Fr> l_apos, r_apos; // - z \cdot 1^D, y_bold ^ D \circuit (z\cdot 1^d) + \sum ...
  std::vector<Fr> sum_z_pow_dj;   // ∑ j in I {z^(1+j) \cdot d_j}
  
  std::vector<Fr> poly_fzy;
  std::vector<Fr> poly_fzy_1, poly_fzy_2;

  base = 0;
  sum_z_pow_dj = std::vector<Fr>(_D, Fr::one());
  for(auto j: this->_ranges.get_indexs()) {
    bit_len = this->_ranges.get_range(j).get_bit_len();
    sum_z_pow_dj[base] = z_powers_0[j + 1];
    for(i = 1; i < bit_len; i++) {
      sum_z_pow_dj[base + i] = sum_z_pow_dj[base + i - 1] * two;
    }
    base += bit_len;
  }
  sum_z_pow_dj = hadmard_product(y_inv_powers_0, sum_z_pow_dj);

  l_apos = std::vector<Fr>(_D, - _z * Fr::one());
  r_apos = std::vector<Fr>(_D, _z * Fr::one());
  r_apos = vector_add(r_apos, sum_z_pow_dj);

  poly_fzy = poly_merge(l_apos, 1, r_apos, _N + 2, _N + _D + 2);

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
  poly_fh = shift_and_reverse(hadmard_product(y_inv_powers_0, ss_inv), _N + 2);

#ifdef CRED_DEBUG
  assert(poly_fg.size() == _D + 1);
  assert(poly_fh.size() == _N + _D + 2);
#endif

  agenda.create_item("Prove::kzg_prove");
  pi._pi_kzg_fzy = kzg_prove(_srs, poly_fzy, poly_fzy.size(), pi._commit_fzy);
  pi._pi_kzg_fg = kzg_prove(_srs, poly_fg, poly_fg.size(), pi._commit_fg);
  pi._pi_kzg_fh = kzg_prove(_srs, poly_fh, poly_fh.size(), pi._commit_fh);
  agenda.mark_item_end("Prove::kzg_prove");

#ifdef CRED_DEBUG
  G1 g, h;
  g = G1::zero();
  h = G1::zero();
  for(i = 0; i < g_vec.size(); i++) {
    g = g + ss[i] * g_vec[i];
    h = h + ss[i].inverse() * h_vec[i];
  }
  assert(pi._commit_fg == g);
  assert(pi._commit_fh == h);

  G1 F;
  h_vec.clear();    // reset h_vec
  for(i = 0; i < this->_D; i++) {
    h_vec.push_back(this->_srs.get_g1_beta_exp(this->_N+2+i));
  }

  F = G1::zero();
  for(size_t i = 0; i < _D; i++) {
    F = F + l_apos[i] * g_vec[i] + r_apos[i] * h_vec[i];
  }
  assert(pi._commit_fzy == F);
#endif

  libff::leave_block("Prover::prove_improved()");

  agenda.mark_item_end("Prove_improved::Round3");
  agenda.mark_item_end("Prove_improved");

  agenda.mark_mem_usage("Prove_improved");
  agenda.mark_proof_size(pi);

  return pi;
}

};