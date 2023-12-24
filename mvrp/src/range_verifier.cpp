#include "structs.hpp"
#include "range_verifier.hpp"

namespace cred {

/************************** Verifier ***************************/
Verifier::Verifier(const CCredSRS& srs, const G1& C, const size_t& set_size, const CRanges& ranges, const IPAProveSystem& ipa_sys) {
  _srs = CCredSRS(srs);
  _set_size = set_size;
  _ranges = CRanges(ranges);
  _ipa_prove_sys = ipa_sys;

  _N = _srs.N;
  _D = _ranges._total_bits;
  _part_sums = _ranges._bits_part_sums;

  _C = G1(C);
}

bool Verifier::verify(const CRangeProof& pi, Agenda& agenda, const bool improved) {
  if(improved)
    return verify_improved(pi, agenda);
  return verify_base(pi, agenda);
}

bool Verifier::verify_base(const CRangeProof& pi, Agenda& agenda) {

  bool result;
  agenda.create_item("Verify_base");

  vector<Fr> transcript_fr;
  vector<G1> transcript_g1;
  vector<G2> transcript_g2;
  vector<GT> transcript_gt;
  vector<Fr> random_challenges;

  Fr _y, _z, _x;  // randoms used in interaction
  
  Fr _delta;
  GT _left_1, _left_2, _left_3, _left;
  GT _right_1, _right_2, _right_3, _right;

#ifdef CRED_DEBUG
  libff::enter_block("Verifier::verify");
#endif

  /** varient of pointproof verify */
  size_t          i;
  Fr              two, two_pow;
  std::vector<Fr> one_vec;    // [1, ..., 1]
  std::vector<Fr> z_powers_0; // index 0 with z^0, [1, z, z^2, ..., z^(D-1)]
  std::vector<Fr> y_powers_0; // index 0 with y^1, [1, y, y^2, ..., y^(D-1)]
  GT beta_N_plus_2_GT;        // (beta ^ (N+2)]_T
  
  two = Fr::one() + Fr::one();
  one_vec = std::vector<Fr>(_D, Fr::one());
  beta_N_plus_2_GT = ReducedPairing(this->_srs.get_g1_beta_exp(_N), this->_srs.get_g2_beta_exp(2));

  const bool pointproof_check = true;
  if(pointproof_check) {

#ifdef CRED_DEBUG
    libff::enter_block("/** varient of pointproof verify */");
#endif

    /** re-compute random y,z,x, copute y_powers, z_powers **/

    transcript_g1.push_back(pi._A);
    transcript_g1.push_back(pi._K);
    random_challenges = gen_random_field_elements_from_transcripts(transcript_g1, transcript_g2, transcript_gt, 2);
    _y = random_challenges[0];
    _z = random_challenges[1];

    transcript_g1.push_back(pi._pi_tilde);
    transcript_gt.push_back(pi._T_1);
    transcript_gt.push_back(pi._T_2);
    _x = gen_random_field_elements_from_transcripts(transcript_g1, transcript_g2, transcript_gt, 1)[0];

    z_powers_0 = vector_powers(_z, _set_size+2); // 1, z, z^2, ... , z^m
    y_powers_0 = vector_powers(_y, _D-1);

    /** compute sum_z_pows_cdot_inner 
     * = ∑j \in I {z ^ (2+j) · <1^nj , [1, 2, ..., 2^(nj-1)]> }
     * = ∑j \in I {z ^ (2+j) · (2^nj - 1) } 
    **/
    Fr sum_z_pows_cdot_inner;  // δ = (z − z^2)·<1^D, y^D> - ∑j \in I {z ^ (2+j) · <1^nj , 2^nj>

    sum_z_pows_cdot_inner = Fr::zero();
    for(size_t j: _ranges.get_indexs()) {
      two_pow = two ^ (_ranges.get_range(j).get_bit_len());
      sum_z_pows_cdot_inner += z_powers_0[j+2] * (two_pow - Fr::one());
    }

    /** compute δ = (z − z^2)·<1^D, y^D> - ∑j \in I {z ^ (2+j) · <1^nj , 2^nj> **/
    _delta = (z_powers_0[1] - z_powers_0[2]) * inner_product(one_vec, y_powers_0) - sum_z_pows_cdot_inner;

    /** compute
     * [left_1]_T = (x · [T1]_T )
     * [left_2]_T = (x^2 · [T1]_T )
     * [left_3]_T = e( C, [∑j in I {z^(1+j) β^(N+1−j)}]_2 )
     * [left]_T = [left_1]_T * [left_2]_T * [left_3]_T
    **/
    G2 sum_z_pows_beta_pows; // sum_z_pows_beta_pows = [∑j in I {z^(1+j) β^(N+1−j)}]_2

    sum_z_pows_beta_pows = G2::zero();
    for(size_t j: _ranges.get_indexs()) {
      sum_z_pows_beta_pows = sum_z_pows_beta_pows + z_powers_0[1+j] * this->_srs.get_g2_beta_exp(_N + 1 - j);
    }

    _left_1 = pi._T_1 ^ _x;
    _left_2 = pi._T_2 ^ (_x ^ 2);
    _left_3 = ReducedPairing(_C, sum_z_pows_beta_pows);
    _left = _left_1 * _left_2 * _left_3;

    /** compute
     * [right_1]_T = [(t_tilde − δ) β ^(N+1)]_T
     * [right_2]_T = [r_x * β ^(N+2)]_T
     * [right_3]_T =  e([pi_tilde]_1, [1]_2)
     * [right]_T = [right_1]_T * [right_2]_T * [right_3]_T
    **/
    _right_1 = _srs.gt ^ (pi._t_tilde - _delta);
    _right_2 = beta_N_plus_2_GT ^ (pi._r_x);
    _right_3 = ReducedPairing(pi._pi_tilde, _srs.get_g2_beta_exp(0));
    // _right_3 = ReducedPairing(pi._pi_tilde, G2::one());
    _right = _right_1 * _right_2 * _right_3;

    result = (_left == _right);

#ifdef CRED_DEBUG
    libff::leave_block("/** varient of pointproof verify */");
#endif

#ifdef CRED_DEBUG
    assert(result);
#endif
  }

  



  /** inner product check for range proof **/
  size_t bit_len, base;
  G1 F, F1, P_apos;
  std::vector<Fr> y_inv_powers_0;   // [1, y, y^2, ..., y^(D-1)]
  std::vector<G1> g_vec;            // [[β^1]_1, ..., [β^D]_1]
  std::vector<G1> h_vec;            // [[y^(1-1) * β^(N+1+1)]_1, ..., [[y^(1-i) * β^(N+1+i)]_1], ..., [y^(1-D) * β^(N+1+D)]_1]
  std::vector<Fr> sum_z_pow_dj;     // ∑ j in I {z^(1+j) \cdot d_j}
  std::vector<Fr> l_apos;           // - z \cdot 1^D
  std::vector<Fr> r_apos;           // [z * y^0, z * y^1, ... , z * y^(D-1)] + sum_z_pow_dj
  
  y_inv_powers_0 = vector_powers(_y.inverse(), _D - 1);

  const bool range_check = true;
  if(range_check) {

#ifdef CRED_DEBUG
    libff::enter_block("/** inner product check for range proof **/");
#endif
    agenda.create_item("/** inner product check for range proof **/");

    for(i = 0; i < this->_D; i++) {
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

      for(i = 1; i < bit_len; i++) {
        sum_z_pow_dj[base + i] = sum_z_pow_dj[base + i - 1] * two;
      }
      base += bit_len;
    }

    F = G1::zero();
    for(i = 0; i < _D; i++) {
      // F = F + l_apos[i] * g_vec[i] + (r_apos[i] + sum_z_pow_dj[i]) * h_vec[i];
      F = F + (r_apos[i] + sum_z_pow_dj[i]) * h_vec[i];
    }
    F = F + (- _z) * sum_up<G1>(g_vec);

    P_apos = pi._A + (_x * pi._K) + F + (-pi._mu) * _srs.get_g1_beta_exp(_N);

    this->_ipa_prove_sys._P = P_apos;
    this->_ipa_prove_sys._c = pi._t_tilde;
    // result &= this->_ipa_prove_sys.IpaVerify(pi._pi_ipa, g_vec, h_vec);
    agenda.create_item("Verify_base::IpaVerify");
    result &= this->_ipa_prove_sys.IpaMultiExpVerify(pi._pi_ipa, g_vec, h_vec); // sometimes wrong
    agenda.mark_item_end("Verify_base::IpaVerify");
    agenda.mark_item_end("/** inner product check for range proof **/");

  
#ifdef CRED_DEBUG
    libff::leave_block("/** inner product check for range proof **/");
#endif

#ifdef CRED_DEBUG
  assert(result);
#endif
  }
  agenda.mark_item_end("Verify_base");
  agenda.mark_mem_usage("Verify_base");

#ifdef CRED_DEBUG
  libff::leave_block("Verifier::verify");
#endif
  return result;
}

bool Verifier::verify_improved(const CRangeProof& pi, Agenda& agenda) {

  bool result;
  agenda.create_item("Verify_improved");

  vector<Fr> transcript_fr;
  vector<G1> transcript_g1;
  vector<G2> transcript_g2;
  vector<GT> transcript_gt;
  vector<Fr> random_challenges;

  Fr _y, _z, _x;  // randoms used in interaction
  
  Fr _delta;
  GT _left_1, _left_2, _left_3, _left;
  GT _right_1, _right_2, _right_3, _right;

#ifdef CRED_DEBUG
  libff::enter_block("Verifier::verify");
#endif

  /** varient of pointproof verify */
  Fr              two, two_pow;
  std::vector<Fr> z_powers_0; // index 0 with z^0, [1, z, z^2, ..., z^(D-1)]
  std::vector<Fr> y_powers_0; // index 0 with y^1, [1, y, y^2, ..., y^(D-1)]
  GT beta_N_plus_2_GT;        // (beta ^ (N+2)]_T
  
  agenda.create_item("beta_N_plus_2_GT");
  two = Fr::one() + Fr::one();
  // beta_N_plus_2_GT = ReducedPairing(this->_srs.get_g1_beta_exp(_N+2), _srs.get_g2_beta_exp(0)));
  beta_N_plus_2_GT = ReducedPairing(this->_srs.get_g1_beta_exp(_N), this->_srs.get_g2_beta_exp(2));
  agenda.mark_item_end("beta_N_plus_2_GT");


  const bool pointproof_check = true;
  if(pointproof_check) {

#ifdef CRED_DEBUG
    libff::enter_block("/** varient of pointproof verify */");
#endif
    /** re-compute random y,z,x, copute y_powers, z_powers **/
    agenda.create_item("Verify_improved::generate_randoms");
    transcript_g1.push_back(pi._A);
    transcript_g1.push_back(pi._K);
    random_challenges = gen_random_field_elements_from_transcripts(transcript_g1, transcript_g2, transcript_gt, 2);
    _y = random_challenges[0];
    _z = random_challenges[1];

    transcript_g1.push_back(pi._pi_tilde);
    transcript_gt.push_back(pi._T_1);
    transcript_gt.push_back(pi._T_2);
    _x = gen_random_field_elements_from_transcripts(transcript_g1, transcript_g2, transcript_gt, 1)[0];

    agenda.mark_item_end("Verify_improved::generate_randoms");

    agenda.create_item("Verify_improved::z_powers_0");
    z_powers_0 = vector_powers(_z, _set_size+2); // 1, z, z^2, ... , z^m
    // y_powers_0 = vector_powers(_y, _D-1);
    agenda.mark_item_end("Verify_improved::z_powers_0");

    /** compute sum_z_pows_cdot_inner 
     * = ∑j \in I {z ^ (2+j) · <1^nj , [1, 2, ..., 2^(nj-1)]> }
     * = ∑j \in I {z ^ (2+j) · (2^nj - 1) } 
    **/
    Fr sum_z_pows_cdot_inner;  // δ = (z − z^2)·<1^D, y^D> - ∑j \in I {z ^ (2+j) · <1^nj , 2^nj>

    agenda.create_item("Verify_improved::sum_z_pows_cdot_inner");
    sum_z_pows_cdot_inner = Fr::zero();
    for(size_t j: _ranges.get_indexs()) {
      two_pow = two ^ (_ranges.get_range(j).get_bit_len());
      sum_z_pows_cdot_inner += z_powers_0[j+2] * (two_pow - Fr::one());
    }
    agenda.mark_item_end("Verify_improved::sum_z_pows_cdot_inner");

    /** compute δ = (z − z^2)·<1^D, [1,...,y^(D-1)]> - ∑j \in I {z ^ (2+j) · <1^nj , 2^nj> **/
    /** <1^D, [1,...,y^(D-1)]> = (y^D - 1) / (y-1) */
    agenda.create_item("Verify_improved::_delta");
    Fr inner_equivalent = _y ^ (_D);
    inner_equivalent   -= Fr::one();
    inner_equivalent   *= (_y - Fr::one()).inverse();
    _delta = (z_powers_0[1] - z_powers_0[2]) * inner_equivalent - sum_z_pows_cdot_inner;
    agenda.mark_item_end("Verify_improved::_delta");

    /** compute
     * [left_1]_T = (x · [T1]_T )
     * [left_2]_T = (x^2 · [T1]_T )
     * [left_3]_T = e( C, [∑j in I {z^(1+j) β^(N+1−j)}]_2 )
     * [left]_T = [left_1]_T * [left_2]_T * [left_3]_T
    **/
    G2 sum_z_pows_beta_pows; // sum_z_pows_beta_pows = [∑j in I {z^(1+j) β^(N+1−j)}]_2
    agenda.create_item("Verify_improved::sum_z_pows_beta_pows");

    sum_z_pows_beta_pows = G2::zero();
    for(size_t j: _ranges.get_indexs()) {
      sum_z_pows_beta_pows = sum_z_pows_beta_pows + z_powers_0[1+j] * this->_srs.get_g2_beta_exp(_N + 1 - j);
    }
    agenda.mark_item_end("Verify_improved::sum_z_pows_beta_pows");

    agenda.create_item("Verify_improved::_left");
    _left_1 = pi._T_1 ^ _x;
    _left_2 = pi._T_2 ^ (_x ^ 2);
    _left_3 = ReducedPairing(_C, sum_z_pows_beta_pows);
    _left = _left_1 * _left_2 * _left_3;
    agenda.mark_item_end("Verify_improved::_left");

    /** compute
     * [right_1]_T = [(t_tilde − δ) β ^(N+1)]_T
     * [right_2]_T = [r_x * β ^(N+2)]_T
     * [right_3]_T =  e([pi_tilde]_1, [1]_2)
     * [right]_T = [right_1]_T * [right_2]_T * [right_3]_T
    **/
    agenda.create_item("Verify_improved::_right");
    _right_1 = _srs.gt ^ (pi._t_tilde - _delta);
    _right_2 = beta_N_plus_2_GT ^ (pi._r_x);
    // _right_3 = ReducedPairing(pi._pi_tilde, G2::one());
    _right_3 = ReducedPairing(pi._pi_tilde, _srs.get_g2_beta_exp(0));
    _right = _right_1 * _right_2 * _right_3;
    agenda.mark_item_end("Verify_improved::_right");

    result = (_left == _right);

#ifdef CRED_DEBUG
    libff::leave_block("/** varient of pointproof verify */");
#endif

#ifdef CRED_DEBUG
  assert(result);
#endif
  }

  
  /** inner product check for range proof **/
  G1 F, g, h, P_apos;
  bool kzg_result_fzy, kzg_result_fg, kzg_result_fh;
  std::vector<G1> g1tuple;
  std::vector<G2> g2tuple;
  CKzgKey kzg_ck;

  const bool range_check = true;
  if(range_check) {

#ifdef CRED_DEBUG
    libff::enter_block("/** inner product check for range proof **/");
#endif
    agenda.create_item("/** inner product check for range proof **/");

    F = pi._commit_fzy;
    // F = pi._commit_fzy_1 + pi._commit_fzy_2;
    P_apos = pi._A + (_x * pi._K) + F + (-pi._mu) * _srs.get_g1_beta_exp(_N);
    this->_ipa_prove_sys._P = P_apos;

    this->_ipa_prove_sys._c = pi._t_tilde;
    
    g = pi._commit_fg;
    h = pi._commit_fh;
    agenda.create_item("Verify_improved::IpaVerify");
    result &= this->_ipa_prove_sys.IpaMultiExpVerify(pi._pi_ipa, g, h);
    agenda.mark_item_end("Verify_improved::IpaVerify");

    g1tuple = {_srs.get_g1_beta_exp(0)};
    g2tuple = {_srs.get_g2_beta_exp(0), _srs.get_g2_beta_exp(1)};
    kzg_ck = CKzgKey(g1tuple, g2tuple);

#ifdef CRED_DEBUG
    libff::enter_block("/** polynomial check **/");
#endif
    agenda.create_item("Verify_improved::kzg_vfyeval");
    kzg_result_fzy = kzg_vfyeval(kzg_ck, pi._commit_fzy, pi._pi_kzg_fzy);
    kzg_result_fg  = kzg_vfyeval(kzg_ck, pi._commit_fg, pi._pi_kzg_fg);
    kzg_result_fh  = kzg_vfyeval(kzg_ck, pi._commit_fh, pi._pi_kzg_fh);
    agenda.mark_item_end("Verify_improved::kzg_vfyeval");

#ifdef CRED_DEBUG
    libff::leave_block("/** polynomial check **/");
#endif

    result &= kzg_result_fzy & kzg_result_fg & kzg_result_fh;
    // result &= kzg_result_fzy;

#ifdef CRED_DEBUG
    libff::leave_block("/** inner product check for range proof **/");
#endif
    agenda.mark_item_end("/** inner product check for range proof **/");

#ifdef CRED_DEBUG
  assert(result);
#endif
  }

  agenda.mark_item_end("Verify_improved");
  agenda.mark_mem_usage("Verify_improved");


#ifdef CRED_DEBUG
  libff::leave_block("Verifier::verify");
#endif

  return result;
}

};
