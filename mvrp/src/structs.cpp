#include <vector>
#include <ostream>
#include <iterator>
#include <sstream>
#include <openssl/sha.h>
#include "structs.hpp"

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
  // this->g1_base_ = G1::G1_one;
  // this->g2_base_ = G2::G2_one;
  // this->gt_base_ = ReducedPairing(G1::G1_one, G2::G2_one);

  this->g1_base_ = G1::random_element();
  this->g2_base_ = G2::random_element();
  this->gt_base_ = ReducedPairing(g1_base_, g2_base_);

  this->N = N;
  this->g1s = std::vector<G1>(2*N - 1);
  this->g2s = std::vector<G2>(N);
  this->beta_powers = vector_powers(beta, 2*N); // index i with exp i

  for(size_t i = 1; i <= N; i++) {
    this->g1s[i-1] = this->beta_powers[i] * g1_base_; // index 0 ~ N-1 with exp 1 ~ N
    this->g2s[i-1] = this->beta_powers[i] * g2_base_; // index 0 ~ N-1 with exp 1 ~ N
  }
  this->gt = gt_base_ ^ this->beta_powers[N+1];  // now exp is beta^(N+1)

  for(size_t i = N+2; i <= 2*N; i++) {
    this->g1s[i-2] = this->beta_powers[i] * g1_base_; // index N ~ 2N - 2 with exp N+2 ~ 2N
  }
}

CredSRS::CredSRS(const CredSRS& other) {
  this->g1_base_ = other.g1_base_;
  this->g2_base_ = other.g2_base_;
  this->gt_base_ = other.gt_base_;
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
    // printf("part_sum[%ld] = %ld\n", i, part_sum);

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

vector<Fr> generate_random_fr_vec(const vector<G1>& vec_1, const vector<G2>& vec_2, const vector<GT>& vec_T, const size_t out_len) {
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

vector<Fr> generate_random_fr_vec(const vector<Fr>& vec_fr, const vector<G1>& vec_1, const vector<G2>& vec_2, const vector<GT>& vec_T, const size_t out_len) {
#ifdef CRED_DEBUG
  assert((vec_fr.size() + vec_1.size()) > 0 && out_len > 0);
#endif

  std::string str;
  vector<Fr>  ret(out_len);
  std::stringstream ss;

  for(auto fr: vec_fr) {
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
// libff::bit_vector generate_hash_bits(string& num_string) {
//   std::vector<size_t> v_int;
//   for (size_t i = 0; i < num_string.size(); i++) {
//     v_int.push_back((size_t)((size_t)num_string.at(i) - 48));
//   }

//   const size_t repack = 1; //what is the use of this repack variable?
//   std::vector<Fr> v_fr = libff::pack_int_vector_into_field_element_vector<Fr>(v_int, repack);
//   libff::bit_vector v_bits_in = libff::convert_field_element_vector_to_bit_vector<Fr>(v_fr);
//   v_bits_in.resize(libsnark::SHA256_block_size);
//   libff::bit_vector v_bits_out = libsnark::sha256_two_to_one_hash_gadget<Fr>::get_hash(v_bits_in);
//   return v_bits_out;
// }

// Fr generate_random(string s) {
//   libff::bit_vector bits;

//   bits = generate_hash_bits(s);
//   bits.resize(Fr::num_bits, true);
//   return libff::convert_bit_vector_to_field_element<Fr>(bits);
// }

// y = hash(A, K), z = hash(y)
void generate_random_y_z(const G1 &A, const G1 &K, Fr &y, Fr &z) {
  
  std::string str;
  libff::bit_vector bits;
  str  = A.coord[2].toString(10);
  str += K.coord[2].toString(10);

  // bits = generate_hash_bits(str);
  // bits.resize(Fr::num_bits, true);
  // y = libff::convert_bit_vector_to_field_element<Fr>(bits);
  y = generate_random_sha256(str);

  // bits.resize(libsnark::SHA256_block_size, true);
  // bits = libsnark::sha256_two_to_one_hash_gadget<Fr>::get_hash(bits);
  // bits.resize(Fr::num_bits, true);
  // z = libff::convert_bit_vector_to_field_element<Fr>(bits);
  std::stringstream ss;
  ss << y.as_bigint();
  z = generate_random_sha256(ss.str());
}

// x = hash(pi_tilde, T1, T2)
void generate_random_x(const G1 &pi_tilde, const GT &T1, const GT &T2, Fr &x) {
  std::string str;
  libff::bit_vector bits;
  str  = pi_tilde.coord[2].toString(10);

  // bits = generate_hash_bits(str);
  // bits.resize(Fr::num_bits, true);
  // x = libff::convert_bit_vector_to_field_element<Fr>(bits);
  x = generate_random_sha256(str);
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

std::vector<Fr> vec_matrix_mult(const vector<Fr>& vec, const vector<vector<Fr>>& matrix, const bool vec_left) {
  size_t rows = matrix.size();
  size_t cols = matrix[0].size();
  std::vector<Fr> result;  // 初始化结果向量为零向量
  if (vec_left) {
    // 向量在矩阵左侧
#ifdef CRED_DEBUG
    assert(vec.size() == matrix.size() && vec.size() > (size_t)0);
#endif
    result.resize(matrix[0].size());
    for (size_t i = 0; i < cols; ++i) {
      for (size_t j = 0; j < rows; ++j) {
        result[i] = result[i] + vec[j] * matrix[j][i];
      }
    }
  } else {
    // 向量在矩阵右侧
#ifdef CRED_DEBUG
    assert(matrix.size() > (size_t)0 && vec.size() == matrix[0].size());
#endif
    result.resize(matrix.size());
    for (size_t i = 0; i < rows; ++i) {
      for (size_t j = 0; j < cols; ++j) {
        result[i] = result[i] + matrix[i][j] * vec[j];
      }
    }
  }

  return result;
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

std::vector<Fr> poly_merge(
  const std::vector<Fr>& poly_1, const size_t base_1, \
  const std::vector<Fr>& poly_2, const size_t base_2, \
  const size_t n ) 
{
  
  size_t n1, n2;
  std::vector<Fr> result;

  n1 = poly_1.size();
  n2 = poly_2.size();
  result = std::vector<Fr>(n+1, Fr::zero());

  for(size_t i = 0; i < n1; i++) {
    result[n - (base_1 + i)] += poly_1[i];
  }
  for(size_t i = 0; i < n2; i++) {
    result[n - (base_2 + i)] += poly_2[i];
  }
  return result;
}

void test_poly_merge() {
  size_t n1, n2, n;
  n1 = (size_t)3;
  n2 = (size_t)5;
  n  = (size_t)9;

  Fr zero = Fr::zero();
  Fr one  = zero + Fr::one();
  Fr two  = one  + Fr::one();
  Fr three  = two + Fr::one();
  Fr four  = three + Fr::one();
  Fr five  = four + Fr::one();
  Fr six  = five + Fr::one();
  Fr seven  = six + Fr::one();
  Fr eight  = seven + Fr::one();
  Fr nine  = eight + Fr::one();

  size_t base_1, base_2;
  base_1 = 1;
  base_2 = 5;
  std::vector<Fr> poly_1 = {one, two, three};
  std::vector<Fr> poly_2 = {five, six, seven, eight, nine};
  std::vector<Fr> poly = {nine, eight, seven, six, five, zero, three, two, one, zero};
  std::vector<Fr> poly_apos = poly_merge(poly_1, base_1, poly_2, base_2, n);

  std::ostream_iterator<Fr>   outFr(cout, ", ");

  libff::print_header("poly");
  copy(poly.begin(), poly.end(), outFr);
  cout << endl;

  libff::print_header("poly_apos");
  copy(poly_apos.begin(), poly_apos.end(), outFr);
  cout << endl;

  assert(vector_all_zero(vector_sub(poly, poly_apos)));
}

bool vector_all_zero(const std::vector<Fr>& v) {
  for(auto item: v) {
    if(item != Fr::zero()) 
      return false;
  }
  return true;
}

std::vector<Fr> multi_exponentiation_nlogn(const std::vector<Fr>& randoms, const std::vector<Fr>& randoms_inv, const size_t recursion_time) {
  Fr x, x_inv;
  size_t step, interval, i, j, k, len;
  std::vector<Fr> ss;

  len = 1 << recursion_time;
  ss = std::vector<Fr>(len, Fr::one());

  for(i = recursion_time; i > 0; i--) {
    step     = (size_t)1 << i;
    interval = step / 2;
    x        = randoms[recursion_time-i];
    x_inv    = randoms_inv[recursion_time-i];
    for(j = 0; j < len; j += step) {
      for(k = j; k < j+interval; k++) {
        ss[k]          *= x_inv;
        ss[k+interval] *= x;
      }
    }
  }
  return ss;
}

std::vector<Fr> multi_exponentiation_n(const std::vector<Fr>& randoms, const size_t recursion_time) {
  size_t i, j, len, base;
  std::vector<Fr> ss;

  len = 1 << recursion_time;
  ss = std::vector<Fr>(len, Fr::one());
  for(i = 0; i < recursion_time; i++) {
    ss[0] *= randoms[i];
  }
  ss[0] = ss[0].inverse();

  for(i = 1; i <= recursion_time; i++) {
    base = 1 << (i-1);
    for(j = 0; j < base; j++) {
      ss[base+j] = ss[j] * (randoms[recursion_time-i] * randoms[recursion_time-i]);
    }
  }
  return ss;
}

AgendaItem::AgendaItem(const string& label) {
  this->label = label;
  this->enter_time      = libff::get_nsec_time();
  this->enter_cpu_time  = libff::get_nsec_cpu_time();
  this->leave_time        = 0;
  this->leave_cpu_time    = 0;

  this->times_enter = 1;
  this->times_leave = 0;
};

void AgendaItem::mark_enter() {
  this->enter_time     += libff::get_nsec_time();
  this->enter_cpu_time += libff::get_nsec_cpu_time();

  this->times_enter += 1;
}

void AgendaItem::mark_leave() {
  this->leave_time     += libff::get_nsec_time();
  this->leave_cpu_time += libff::get_nsec_cpu_time();

  this->times_leave += 1;
}

void AgendaItem::print(long long& start_time, long long& start_cpu_time) {
  printf("%-35s\t", this->label.c_str());
  assert(this->times_enter == this->times_leave);
  
  // long long time_from_start = leave_time - start_time;
  long long time_from_last = (leave_time - enter_time) / this->times_enter;

  // long long cpu_time_from_start = leave_cpu_time - start_cpu_time;
  long long cpu_time_from_last = (leave_cpu_time - enter_cpu_time) / this->times_enter;

  if (time_from_last != 0) {
      double parallelism_from_last = 1.0 * cpu_time_from_last / time_from_last;
      printf("[%0.8fs x%0.2f]\n", time_from_last * 1e-9, parallelism_from_last);
  } else {
      printf("[             ]\n");
  }
}

Agenda::Agenda(const string mode, const size_t N, const size_t D, const size_t num, const size_t bit_len) {
  start_time = libff::get_nsec_time();
  start_cpu_time = libff::get_nsec_cpu_time();

  this->mode = mode;
  this->N = N;
  this->D = D;
  this->num = num;
  this->bit_len = bit_len;
  items.clear();
};

void Agenda::create_item(const string& label) {
  if(this->items.find(label) == this->items.end()) {
    label_records.push_back(label);

    // long long now = libff::get_nsec_time();
    AgendaItem item = AgendaItem(label);
    items[label] = item;
  } else {
    items[label].mark_enter();
  }


}

void Agenda::mark_item_end(const string& label) {
  assert(this->items.find(label) != this->items.end());

  items[label].mark_leave();
}

string Agenda::to_string_rough() {
  char ret[256];
  if(mode == "base") {
    sprintf(ret, "%2lu bit \\times\\thinspace %4lu & %4lu & %0.1f & %0.1f\n", bit_len, num, proof_size, items["Prove_base"].duration() * 1000, items["Verify_base"].duration() * 1000);
  } else {
    sprintf(ret, "%2lu bit \\times\\thinspace %4lu & %4lu & %0.1f & %0.1f\n", bit_len, num, proof_size, items["Prove_improved"].duration() * 1000, items["Verify_improved"].duration() * 1000);
  }
  return string(ret);
}

string Agenda::ipa_record() {
  char ret[256];
  assert(mode == "base");
  sprintf(ret, "%2lu bit \\times\\thinspace %4lu & %4lu & %0.1f & %0.1f\n", bit_len, num, ipa_size, items["Prove::IpaProve"].duration() * 1000, items["Verify_base::IpaVerify"].duration() * 1000);
  return string(ret);
}

string Agenda::bullet_record() {
  char ret[256];
  assert(mode == "base");
  float bullet_prove_t = 1000*(items["bullet::Round1"].duration() + items["bullet::Round2_plus_pi_tilde"].duration() \
                               - items["Prove::_pi_tilde"].duration() + items["bullet::Round3"].duration());
  sprintf(ret, "%2lu bit \\times\\thinspace %4lu & %4lu & %0.1f & %0.1f\n", bit_len, num, ipa_size, bullet_prove_t, items["Verify_base::IpaVerify"].duration() * 1000);
  return string(ret);
}

void Agenda::print() {
  printf("Mode: %s\n", mode.c_str());
  printf("Args: N_%lu, D_%lu, num_%lu, bitlen_%lu\n", N, D, num, bit_len);
  printf("ProofSize: %lu\n", proof_size);
  for(auto label : label_records) {
    items[label].print(this->start_time, this->start_cpu_time);
  }

  printf("\n\n");
  for(auto item: v_sizes) {
    printf("* Peak vsize (physical memory+swap) in mebibytes (%s): %lu\n", item.first.c_str(), item.second >> 20);
  }
}

void Agenda::write_file(const char* filename) {
  FILE* fp = freopen(filename,"w", stdout);
  if(fp == NULL) cout<<"error"<<endl;

  this->print();
  fclose(fp);

  fp = freopen("/dev/tty","w",stdout); // 恢复默认输出
  // fclose(fp);
}

bool write_file(const string& file_path, const string& content) {
  FILE* fp = fopen(file_path.c_str(),"w");
  if(fp == NULL) cout<<"error"<<endl;

  // printf("sizeof(content) %lu\n", sizeof(char) * content.size());
  fwrite(content.c_str(), sizeof(char) * content.size(), 1, fp);
  fclose(fp);

  return true;
}


};
