#ifndef CRED_STRUCTS_HPP
#define CRED_STRUCTS_HPP
#include <vector>
#include <map>
#include <tuple>
#include "libff/common/profiling.hpp"
#include "libff/algebra/curves/public_params.hpp"
#include "libff/common/default_types/ec_pp.hpp"
#include <xassert/XAssert.h>


using namespace std;

namespace cred {

// Type of group G1
using G1 = typename libff::default_ec_pp::G1_type;
// Type of group G2
using G2 = typename libff::default_ec_pp::G2_type;
// Type of group GT (recall pairing e : G1 x G2 -> GT)
using GT = typename libff::default_ec_pp::GT_type;
// Type of the finite field "in the exponent" of the EC group elements
using Fr = typename libff::default_ec_pp::Fp_type;

// Pairing function, takes an element in G1, another in G2 and returns the one in GT
//using libff::default_ec_pp::reduced_pairing;
//using ECPP::reduced_pairing;
// Found this solution on SO: https://stackoverflow.com/questions/9864125/c11-how-to-alias-a-function
// ('using ECPP::reduced_pairing;' doesn't work, even if you expand ECPP)
template<typename ... Args>
auto ReducedPairing(
        Args&&... args) -> decltype(libff::default_ec_pp::reduced_pairing(std::forward<Args>(args)...)) {
    return libff::default_ec_pp::reduced_pairing(std::forward<Args>(args)...);
}

class CredSRS
{
public:
  std::vector<G1> g1s;
  std::vector<G2> g2s;
  std::vector<Fr> beta_powers;
  GT gt;
  size_t N;

  CredSRS() {};
  CredSRS(size_t N);
  CredSRS(const CredSRS& other);

  const G1& get_g1_beta_exp(size_t);
  const G2& get_g2_beta_exp(size_t);

  static G1& get_g1_beta_exp(CredSRS& self, size_t idx);
  static G2& get_g2_beta_exp(CredSRS& self, size_t idx);
};

class Range {
public:
  Fr low, high;   // 下限, 上限
  size_t bit_len; // 上限的比特位数

  Range() {};
  Range(Fr low, Fr high, size_t bit_len);
  Range(const Range& other);

  inline const size_t get_bit_len() const {
    return this->bit_len;
  };
};

class Ranges {
public:
  std::vector<size_t> _indexs;
  std::map<size_t, Range> _range_map;

  size_t _total_bits;
  std::map<size_t, size_t> _bits_part_sums;

  Ranges() {};
  Ranges(const std::vector<size_t>& indexs, const std::vector<Range>& range_vec);
  Ranges(const Ranges& other);

  inline const std::vector<size_t>& get_indexs() {
    return this->_indexs;
  };
  
  inline const Range& get_range(size_t i) {
    return this->_range_map[i];
  };

  inline const size_t get_part_sum(size_t i) {
    return this->_bits_part_sums[i];
  };

private:
  size_t total_bits();
  std::map<size_t, size_t> bits_part_sums();
};

class SET
{
public:
  size_t m; // set size m
  std::vector<Fr> set_values;

  SET() {};
  SET(std::vector<Fr> set_values);
  SET(const SET& other);
  const Fr& get_set_value(size_t idx);
};

class IPAProof {
public:
  std::vector<G1> L_vec;
  std::vector<G1> R_vec;
  Fr a, b;
};

class IPAProofOneRecursion {
public:
  G1 L;
  G1 R;
  std::vector<Fr> a_vec;
  std::vector<Fr> b_vec;
};

class CredProof {
public:
  // round 1
  G1 _A, _K;

  // round 2
  G1 _F;
  GT _T_1, _T_2;

  // round 3
  Fr _mu;   // r_1 + r_2 * x
  Fr _r_x;  // r_3 * x + r_4 * x^2
  Fr _t_tilde;  // t(x)
  Fr _f_s;      // q_s = f(s)
  G1 _pi_tilde; // product_(j \in _I) pi_j
  G1 _pi_Q;     // [q(beta) * beta^(N+2-D)]_1

  // G1 _P_to_del; // this should be computed by the verifier himself, and should not be part of proof. 
                // But in the code, the verifier can't compute the right P now, so I passed the prover
                // computed P.
  IPAProof _pi_ipa;   // 
};

// 使用位运算判断是否只有一个位为1
bool isPowerOfTwo(size_t n);

void bits_padding(libff::bit_vector& v);

Fr generate_random(string s);

void generate_random_y_z(const G1 &A, const G1 &K, Fr &y, Fr &z);

void generate_random_x(const G1 &pi_tilde, const GT &T1, const GT &T2, Fr &x);

Fr inner_product(const std::vector<Fr> &v1, const std::vector<Fr> &v2);

// need to check if the relation of v1, v2, len
Fr inner_product(std::vector<Fr>::const_iterator v1_iter, std::vector<Fr>::const_iterator v2_iter, size_t len);

std::vector<Fr> hadmard_product(const std::vector<Fr> &v1, const std::vector<Fr> &v2);

std::vector<Fr> numerical_mult(const Fr &factor, const std::vector<Fr> &v);

std::vector<Fr> numerical_mult(const Fr &factor, std::vector<Fr>::const_iterator v_begin, size_t size);

std::vector<Fr> vector_add(const std::vector<Fr> &v1, const std::vector<Fr> &v2);

std::vector<Fr> vector_sub(const std::vector<Fr> &v1, const std::vector<Fr> &v2);

std::vector<Fr> vector_neg(const std::vector<Fr> &v);

std::vector<Fr> vector_powers(const Fr &x, size_t max_exp, bool zero_exp = true);
};
#endif