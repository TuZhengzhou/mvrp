#include <openssl/sha.h>
#include <iostream>
#include "basic_types.hpp"
#include "math_operations.hpp"

using std::cout;
using std::endl;

namespace cred {

CCredSRS::CCredSRS(size_t N) {
  if(N == 0) 
    return;

  Fr beta = Fr::random_element();
  G1 g1_base;
  G2 g2_base;
  GT gt_base;

  g1_base = G1::random_element();
  g2_base = G2::random_element();
  gt_base = ReducedPairing(g1_base, g2_base);

  this->N = N;
  this->g1s = std::vector<G1>(2*N + 1);
  this->g2s = std::vector<G2>(N + 1);
  this->beta_powers = vector_powers(beta, 2*N); // index i with exp i

  for(size_t i = 0; i <= N; i++) {
    this->g1s[i] = this->beta_powers[i] * g1_base; // index 0 ~ N-1 with exp 1 ~ N
    this->g2s[i] = this->beta_powers[i] * g2_base; // index 0 ~ N-1 with exp 1 ~ N
  }
  // g1s[N+1] will be used since kzg_prove has no index check, so we set it to zero
  this->g1s[N+1] = G1::zero();
  this->gt = gt_base ^ this->beta_powers[N+1];  // now exp is beta^(N+1)
  // tmpararily set index N+1 with exp N+1 to test the correctness of the code
  // this->g1s[N+1] = this->beta_powers[N+1] * g1_base; // index N+1 with exp N+1

  for(size_t i = N+2; i <= 2*N; i++) {
    this->g1s[i] = this->beta_powers[i] * g1_base; // index N ~ 2N - 2 with exp N+2 ~ 2N
  }
}

CCredSRS::CCredSRS(const CCredSRS& other) {
  this->g1s = std::vector<G1>(other.g1s);
  this->g2s = std::vector<G2>(other.g2s);
  this->beta_powers = std::vector<Fr>(other.beta_powers);
  this->gt = other.gt;
  this->N = other.N;
}

const G1& CCredSRS::get_g1_beta_exp(size_t exp) const {
  assert(exp != N+1 && exp <= 2*N);
  return this->g1s[exp];
}

const G2& CCredSRS::get_g2_beta_exp(size_t exp) const {
  assert(exp <= N);
  return this->g2s[exp];
}

/************************** CRange ***************************/
CRange::CRange(Fr low, Fr high, size_t bit_len) {
  this->low = low;
  this->high = high;
  this->bit_len = bit_len;
}

CRange::CRange(const CRange& other) {
  this->low = other.low;
  this->high = other.high;
  this->bit_len = other.bit_len;
}

size_t CRanges::total_bits() {
  size_t ans = 0;
  for(const std::pair<size_t, CRange>& item: this->_range_map) {
    ans = ans + item.second.get_bit_len();
  }
  return ans;
}

std::map<size_t, size_t> CRanges::bits_part_sums() {
  std::map<size_t, size_t> ans;
  if(this->_indexs.size() == 0)
    return ans;

  size_t part_sum = 0;
  for(size_t i: this->_indexs) {
    ans.insert(std::pair<size_t, size_t>(i, part_sum));
    // printf("part_sum[%ld] = %ld\n", i, part_sum);

    part_sum += this->_range_map[i].get_bit_len();
  }
  return ans;
}

CRanges::CRanges(const std::vector<size_t>& indexs, const std::vector<CRange>& range_vec) {

  assert(indexs.size() == range_vec.size());

  size_t len = indexs.size();
  this->_indexs = indexs;

  for(size_t i = 0; i < len; i++) {
    this->_range_map[indexs[i]] = range_vec[i];
  }
  this->_total_bits = this->total_bits();
  this->_bits_part_sums = this->bits_part_sums();
}

CRanges::CRanges(const CRanges& other) {
  this->_indexs = other._indexs;
  this->_range_map = other._range_map;
  this->_total_bits = other._total_bits;
  this->_bits_part_sums = other._bits_part_sums;
}

/************************** CSet ***************************/

CSet::CSet(std::vector<Fr> set_values) {
  this->m = set_values.size();
  this->set_values = std::vector<Fr>(set_values);
}

CSet::CSet(const CSet& other) {
  this->m = other.m;
  this->set_values = std::vector<Fr>(other.set_values);
}

const Fr& CSet::get_set_value(size_t idx) {
  assert(idx >= 1 && idx <= this->m);
  return this->set_values[idx-1];  // index 0 ~ N-1 with exp 1 ~ N
}

void toString(std::vector<Fr>&v) {
  std::stringstream ss;
  for(auto item: v) {
    ss << item.as_ulong() << " ";
  }
  printf("%s\n", ss.str().c_str());
}

/************************** utils ***************************/
bool is_power_of_two(size_t n) {
  // 零不是二的幂次
  if (n == (size_t)0) {
      return false;
  }
  return (n & (n - 1)) == 0;
}

std::string hexToChar(const char c) {
    switch(tolower(c))
    {
        case '0': return "0000";
        case '1': return "0001";
        case '2': return "0010";
        case '3': return "0011";
        case '4': return "0100";
        case '5': return "0101";
        case '6': return "0110";
        case '7': return "0111";
        case '8': return "1000";
        case '9': return "1001";
        case 'a': return "1010";
        case 'b': return "1011";
        case 'c': return "1100";
        case 'd': return "1101";
        case 'e': return "1110";
        case 'f': return "1111";
    }
    return "";
}

libff::bit_vector hexToBin(const std::string& str) {
    libff::bit_vector res(str.size() * 4);
    for(size_t i = 0; i < str.size(); i++) {
        std::string hexItem = hexToChar(str[i]);
        res[i*4+0] = hexItem[0] == '1' ? true : false;
        res[i*4+1] = hexItem[1] == '1' ? true : false;
        res[i*4+2] = hexItem[2] == '1' ? true : false;
        res[i*4+3] = hexItem[3] == '1' ? true : false;
    }
    return res;
}

Fr generate_random_sha256(const std::string& pre_image) {

  unsigned char hash[SHA256_DIGEST_LENGTH];
  SHA256((const unsigned char*)pre_image.c_str(), pre_image.length(), hash);

  // 将哈希值打印为十六进制字符串
  std::string hash_hex;
  char buf[3];
  for (int i = 0; i < SHA256_DIGEST_LENGTH; i++) {
      sprintf(buf, "%02x", hash[i]);
      hash_hex += buf;
  }

  libff::bit_vector hash_bits = hexToBin(hash_hex);
  hash_bits.resize(Fr::num_bits, true);
  return libff::convert_bit_vector_to_field_element<Fr>(hash_bits);
}

vector<Fr> gen_random_field_elements_from_transcripts(const vector<G1>& vec_1, const vector<G2>& vec_2, const vector<GT>& vec_T, const size_t out_len) {
#ifdef CRED_DEBUG
  assert(vec_1.size() > 0 && out_len > 0);
#endif

  std::string str;
  vector<Fr>  ret(out_len);
  std::stringstream ss;

  for(auto g1: vec_1) {
    str += g1.coord[2].toString(10);
  }

  ret[0] = generate_random_sha256(str);
  for(size_t i = 1; i < out_len; i++) {
    ss << ret[i-1].as_bigint();
    ret[i] = generate_random_sha256(ss.str());
  }

  return ret;
}

vector<Fr> gen_random_field_elements_from_transcripts(const vector<Fr>& random_challenges, const vector<G1>& vec_1, const vector<G2>& vec_2, const vector<GT>& vec_T, const size_t out_len) {
#ifdef CRED_DEBUG
  assert((random_challenges.size() + vec_1.size()) > 0 && out_len > 0);
#endif

  std::string str;
  vector<Fr>  ret(out_len);
  std::stringstream ss;

  for(auto fr: random_challenges) {
    str += std::to_string(fr.as_ulong());
  }

  for(auto g1: vec_1) {
    str += g1.coord[2].toString(10);
  }

  ret[0] = generate_random_sha256(str);
  for(size_t i = 1; i < out_len; i++) {
    ss << ret[i-1].as_bigint();
    ret[i] = generate_random_sha256(ss.str());
  }

  return ret;
}

bool write_file(const string& file_path, const string& content) {
  FILE* fp = fopen(file_path.c_str(),"w");
  if(fp == NULL) cout<<"error"<<endl;

  // printf("sizeof(content) %lu\n", sizeof(char) * content.size());
  fwrite(content.c_str(), sizeof(char) * content.size(), 1, fp);
  fclose(fp);

  return true;
}

} // namespace cred