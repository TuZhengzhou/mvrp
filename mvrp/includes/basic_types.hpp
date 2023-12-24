#ifndef BASIC_TYPES_HPP
#define BASIC_TYPES_HPP

#include <vector>
#include <string>
#include "libff/common/profiling.hpp"
#include "libff/algebra/curves/public_params.hpp"
#include "libff/common/default_types/ec_pp.hpp"

using std::vector;
using std::string;

namespace cred {

using G1 = typename libff::default_ec_pp::G1_type;
// Type of group G2
using G2 = typename libff::default_ec_pp::G2_type;
// Type of group GT (recall pairing e : G1 x G2 -> GT)
using GT = typename libff::default_ec_pp::GT_type;
// Type of the finite field "in the exponent" of the EC group elements
using Fr = typename libff::default_ec_pp::Fp_type;

class CCredSRS
{
public:
  std::vector<G1> g1s;
  std::vector<G2> g2s;
  std::vector<Fr> beta_powers;
  GT gt;
  size_t N;

  CCredSRS() {};
  CCredSRS(size_t N);
  CCredSRS(const CCredSRS& other);

  const G1& get_g1_beta_exp(size_t) const;
  const G2& get_g2_beta_exp(size_t) const;
};

class CRange {
public:
  Fr low, high;   // 下限, 上限
  size_t bit_len; // 上限的比特位数

  CRange() {};
  CRange(Fr low, Fr high, size_t bit_len);
  CRange(const CRange& other);

  inline const size_t get_bit_len() const { return this->bit_len; };
};

class CRanges {
public:
  std::vector<size_t> _indexs;
  std::map<size_t, CRange> _range_map;

  size_t _total_bits;
  std::map<size_t, size_t> _bits_part_sums;

  CRanges() {};
  CRanges(const std::vector<size_t>& indexs, const std::vector<CRange>& range_vec);
  CRanges(const CRanges& other);

  inline const std::vector<size_t>& get_indexs() { return this->_indexs; };
  
  inline const CRange& get_range(size_t i) { return this->_range_map[i]; };

  inline const size_t get_part_sum(size_t i) { return this->_bits_part_sums[i]; };

private:
  size_t total_bits();
  std::map<size_t, size_t> bits_part_sums();
};

class CSet
{
public:
  size_t m; // set size m
  std::vector<Fr> set_values;

  CSet() {};
  CSet(std::vector<Fr> set_values);
  CSet(const CSet& other);
  const Fr& get_set_value(size_t idx);
};

template<typename ... Args>
auto ReducedPairing(
        Args&&... args) -> decltype(libff::default_ec_pp::reduced_pairing(std::forward<Args>(args)...)) {
    return libff::default_ec_pp::reduced_pairing(std::forward<Args>(args)...);
}

void toString(std::vector<Fr>&v);

// 使用位运算判断是否只有一个位为1
bool is_power_of_two(size_t n);

Fr generate_random_sha256(const std::string& pre_image);

// Fr generate_random(string s);
vector<Fr> gen_random_field_elements_from_transcripts(const vector<G1>& vec_1, const vector<G2>& vec_2, const vector<GT>& vec_T, const size_t out_len);

vector<Fr> gen_random_field_elements_from_transcripts(const vector<Fr>& random_challenges, const vector<G1>& vec_1, const vector<G2>& vec_2, const vector<GT>& vec_T, const size_t out_len);

bool write_file(const string& file_path, const string& content);

}

#endif