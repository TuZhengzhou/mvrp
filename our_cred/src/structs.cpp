#include <vector>
#include <ostream>
#include "structs.hpp"
#include "libsnark/gadgetlib1/gadgets/hashes/sha256/sha256_gadget.hpp"
#include "libsnark/gadgetlib1/gadgets/hashes/sha256/sha256_components.hpp"

using cred::G1;
using cred::G2;
using cred::GT;
using cred::Fr;

/************************** CredSRS ***************************/
namespace cred {

CredSRS::CredSRS(size_t N) {
  if(N == 0) 
    return;

  Fr beta = Fr::random_element();
  G1 g1_one = G1::G1_one;
  G2 g2_one = G2::G2_one;
  GT gt_one = ReducedPairing(g1_one, g2_one);

  this->N = N;
  this->g1s = std::vector<G1>(2*N - 1);
  this->g2s = std::vector<G2>(N);
  this->beta_powers = vector_powers(beta, 2*N); // index i with exp i

  for(size_t i = 1; i <= N; i++) {
    this->g1s[i-1] = this->beta_powers[i] * g1_one; // index 0 ~ N-1 with exp 1 ~ N
    this->g2s[i-1] = this->beta_powers[i] * g2_one; // index 0 ~ N-1 with exp 1 ~ N
  }
  this->gt = gt_one ^ this->beta_powers[N+1];  // now exp is beta^(N+1)

  for(size_t i = N+2; i <= 2*N; i++) {
    this->g1s[i-2] = this->beta_powers[i] * g1_one; // index N ~ 2N - 2 with exp N+2 ~ 2N
  }
}

CredSRS::CredSRS(const CredSRS& other) {
  this->g1s = std::vector<G1>(other.g1s);
  this->g2s = std::vector<G2>(other.g2s);
  this->beta_powers = std::vector<Fr>(other.beta_powers);
  this->gt = other.gt;
  this->N = other.N;
}

const G1& CredSRS::get_g1_beta_exp(size_t exp) {
  assert((exp >= 1 && exp <= N) || (exp >= N+2 && exp <= 2*N));
  if(exp >= 1 && exp <= N) {
    return this->g1s[exp - 1];  // index 0 ~ N-1 with exp 1 ~ N
  } 
  return this->g1s[exp - 2];
}

const G2& CredSRS::get_g2_beta_exp(size_t exp) {
  assert(exp >= 1 && exp <= N);
  return this->g2s[exp - 1];  // index 0 ~ N-1 with exp 1 ~ N
}

/************************** Range ***************************/
Range::Range(Fr low, Fr high, size_t bit_len) {
  this->low = low;
  this->high = high;
  this->bit_len = bit_len;
}

Range::Range(const Range& other) {
  this->low = other.low;
  this->high = other.high;
  this->bit_len = other.bit_len;
}

size_t Ranges::total_bits() {
  size_t ans = 0;
  for(const std::pair<size_t, Range>& item: this->_range_map) {
    ans = ans + item.second.get_bit_len();
  }
  return ans;
}

std::map<size_t, size_t> Ranges::bits_part_sums() {
  std::map<size_t, size_t> ans;
  if(this->_indexs.size() == 0)
    return ans;

  size_t part_sum = 0;
  for(size_t i: this->_indexs) {
    ans.insert(std::pair<size_t, size_t>(i, part_sum));
    printf("part_sum[%ld] = %ld\n", i, part_sum);

    part_sum += this->_range_map[i].get_bit_len();
  }
  return ans;
}

Ranges::Ranges(const std::vector<size_t>& indexs, const std::vector<Range>& range_vec) {

  assert(indexs.size() == range_vec.size());

  size_t len = indexs.size();
  this->_indexs = indexs;

  for(size_t i = 0; i < len; i++) {
    this->_range_map[indexs[i]] = range_vec[i];
  }
  this->_total_bits = this->total_bits();
  this->_bits_part_sums = this->bits_part_sums();
}

Ranges::Ranges(const Ranges& other) {
  this->_indexs = other._indexs;
  this->_range_map = other._range_map;
  this->_total_bits = other._total_bits;
  this->_bits_part_sums = other._bits_part_sums;
}

/************************** SET ***************************/

SET::SET(std::vector<Fr> set_values) {
  this->m = set_values.size();
  this->set_values = std::vector<Fr>(set_values);
}

SET::SET(const SET& other) {
  this->m = other.m;
  this->set_values = std::vector<Fr>(other.set_values);
}

const Fr& SET::get_set_value(size_t idx) {
  assert(idx >= 1 && idx <= this->m);
  return this->set_values[idx-1];  // index 0 ~ N-1 with exp 1 ~ N
}


/************************** utils ***************************/
bool isPowerOfTwo(size_t n) {
  // 零不是二的幂次
  if (n == (size_t)0) {
      return false;
  }
  return (n & (n - 1)) == 0;
}

void bits_padding(libff::bit_vector& v) {
  size_t remainder = v.size() % libsnark::SHA256_block_size;
  v.resize(v.size() + libsnark::SHA256_block_size - remainder, true);
}

libff::bit_vector generate_hash_bits(string& num_string) {
  std::vector<size_t> v_int;
  for (size_t i = 0; i < num_string.size(); i++) {
    v_int.push_back((size_t)((size_t)num_string.at(i) - 48));
  }

  const size_t repack = 1; //what is the use of this repack variable?
  std::vector<Fr> v_fr = libff::pack_int_vector_into_field_element_vector<Fr>(v_int, repack);
  libff::bit_vector v_bits_in = libff::convert_field_element_vector_to_bit_vector<Fr>(v_fr);
  v_bits_in.resize(libsnark::SHA256_block_size);
  libff::bit_vector v_bits_out = libsnark::sha256_two_to_one_hash_gadget<Fr>::get_hash(v_bits_in);
  return v_bits_out;
}

Fr generate_random(string s) {
  libff::bit_vector bits;

  bits = generate_hash_bits(s);
  bits.resize(Fr::num_bits, true);
  return libff::convert_bit_vector_to_field_element<Fr>(bits);
}

// y = hash(A, K), z = hash(y)
void generate_random_y_z(const G1 &A, const G1 &K, Fr &y, Fr &z) {
  
  std::string str;
  libff::bit_vector bits;
  str  = A.coord[2].toString(10);
  str += K.coord[2].toString(10);

  bits = generate_hash_bits(str);
  bits.resize(Fr::num_bits, true);
  y = libff::convert_bit_vector_to_field_element<Fr>(bits);

  bits.resize(libsnark::SHA256_block_size, true);
  bits = libsnark::sha256_two_to_one_hash_gadget<Fr>::get_hash(bits);
  bits.resize(Fr::num_bits, true);
  z = libff::convert_bit_vector_to_field_element<Fr>(bits);
}

// x = hash(pi_tilde, T1, T2)
void generate_random_x(const G1 &pi_tilde, const GT &T1, const GT &T2, Fr &x) {
  std::string str;
  libff::bit_vector bits;
  str  = pi_tilde.coord[2].toString(10);

  bits = generate_hash_bits(str);
  bits.resize(Fr::num_bits, true);
  x = libff::convert_bit_vector_to_field_element<Fr>(bits);
}

Fr inner_product(const std::vector<Fr> &v1, const std::vector<Fr> &v2) {
  assert(v1.size() == v2.size());
  Fr ret = Fr::zero();
  size_t n = v1.size();
  for(size_t i = 0; i < n; i++) {
    ret += v1[i] * v2[i];
  }
  return ret;
}

Fr inner_product(std::vector<Fr>::const_iterator v1_iter, std::vector<Fr>::const_iterator v2_iter, size_t len) {
  Fr ret = Fr::zero();
  for(size_t i = 0; i < len; i++) {
    ret += (*(v1_iter+i)) * (*(v2_iter+i));
  }
  return ret;
}

std::vector<Fr> hadmard_product(const std::vector<Fr> &v1, const std::vector<Fr> &v2) {
  assert(v1.size() == v2.size());
  size_t n = v1.size();
  std::vector<Fr> ret(n);
  for(size_t i = 0; i < n; i++) {
    ret[i] = v1[i] * v2[i];
  }
  return ret;
}

std::vector<Fr> hadmard_product(std::vector<Fr>::const_iterator v1_iter, std::vector<Fr>::const_iterator v2_iter, size_t size) {
  std::vector<Fr> ret(size);
  for(size_t i = 0; i < size; i++) {
    ret[i] = (*(v1_iter+i)) * (*(v2_iter+i));
  }
  return ret;
}


std::vector<Fr> numerical_mult(const Fr &factor, const std::vector<Fr> &v) {
  size_t n = v.size();
  std::vector<Fr> ret(n);
  for(size_t i = 0; i < n; i++) {
    ret[i] = factor * v[i];
  }
  return ret;
}

std::vector<Fr> numerical_mult(const Fr &factor, std::vector<Fr>::const_iterator v_begin, size_t size) {
  std::vector<Fr> ret(size);
  for(size_t i = 0; i < size; i++) {
    ret[i] = factor * (*(v_begin + i));
  }
  return ret;
}


std::vector<Fr> vector_add(const std::vector<Fr> &v1, const std::vector<Fr> &v2) {
  assert(v1.size() == v2.size());
  size_t n = v1.size();
  std::vector<Fr> ret(n);
  for(size_t i = 0; i < n; i++) {
    ret[i] = v1[i] + v2[i];
  }
  return ret;
}

std::vector<Fr> vector_sub(const std::vector<Fr> &v1, const std::vector<Fr> &v2) {
  assert(v1.size() == v2.size());
  size_t n = v1.size();
  std::vector<Fr> ret(n);
  for(size_t i = 0; i < n; i++) {
    ret[i] = v1[i] - v2[i];
  }
  return ret;
}


std::vector<Fr> vector_neg(const std::vector<Fr> &v) {
  Fr zero = Fr::zero();
  std::vector<Fr> ret(v.size());
  for(size_t i = 0; i < v.size(); i++) {
    ret[i] = zero - v[i];
  }
  return ret;
}

std::vector<Fr> vector_powers(const Fr &x, size_t max_exp, bool zero_exp) {

  size_t len = max_exp + (size_t)zero_exp;
  std::vector<Fr> res(len);
  
  res[0] = zero_exp ? Fr::one() : x;
  for (size_t i = 1; i < len; i++) {
    res[i] = res[i-1] * x;
  }

  return res;
}

};
