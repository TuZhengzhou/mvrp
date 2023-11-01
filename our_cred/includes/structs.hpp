#ifndef CRED_STRUCTS_HPP
#define CRED_STRUCTS_HPP
#include <vector>
#include <map>
#include <tuple>
#include "libff/common/profiling.hpp"
#include "libff/algebra/curves/public_params.hpp"
#include "libff/common/default_types/ec_pp.hpp"
#include <xassert/XAssert.h>

#ifndef NO_PROCPS
#include <proc/readproc.h>
#endif

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
  G1 g1_base_;
  G2 g2_base_;
  GT gt_base_;
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

  inline const size_t size() const {
    return 2 * sizeof(Fr)+ 2 * L_vec.size() * sizeof(G1);
  }
};

class IPAProofOneRecursion {
public:
  G1 L;
  G1 R;
  std::vector<Fr> a_vec;
  std::vector<Fr> b_vec;
};

/******************************** Commitment key ********************************/
class kzgCommitKey;
std::ostream& operator<<(std::ostream &out, const kzgCommitKey &ck);
std::istream& operator>>(std::istream &in, kzgCommitKey &ck);

class kzgCommitKey
{
    public:
    std::vector<G1> g1;
    std::vector<G2> g2;

    kzgCommitKey() = default;
    kzgCommitKey& operator=(const kzgCommitKey &other) = default;
    kzgCommitKey(const kzgCommitKey &other) = default;
    kzgCommitKey(kzgCommitKey &&other) = default;
    kzgCommitKey(
        std::vector<G1> &&g1,
        std::vector<G2> &&g2) :
    g1(std::move(g1)),
    g2(std::move(g2))
    {};
    kzgCommitKey(
        std::vector<G1> &g1,
        std::vector<G2> &g2) :
    g1(g1),
    g2(g2)
    {};

    size_t G1_size() const
    {
        return g1.size();
    }

    size_t G2_size() const
    {
        return g2.size();
    }

    size_t GT_size() const
    {
        return 1;
    }

    size_t size_in_bits() const
    {
      size_t ret;
      ret += g1.size() > 0 ? g1.size() * sizeof(g1[0]) : 0;
      ret += g2.size() > 0 ? g2.size() * sizeof(g2[0]) : 0;
      return ret;
    }

    void print_size() const
    {
        libff::print_indent(); printf("* G1 elements in CommitKey: %zu\n", this->G1_size());
        libff::print_indent(); printf("* G2 elements in CommitKey: %zu\n", this->G2_size());
        libff::print_indent(); printf("* Commit Key size in bits: %zu\n", this->size_in_bits());
    }

    bool operator==(const kzgCommitKey &other) const;
    friend std::ostream& operator<< (std::ostream &out, const kzgCommitKey &ck);
    friend std::istream& operator>> (std::istream &in, kzgCommitKey &ck);
};

/******************************** Witness ********************************/

class KZGProof;
std::ostream& operator<< (std::ostream &out, const KZGProof &wit);
std::istream& operator>> (std::istream &in, KZGProof &wit);

class KZGProof
{
    public:
    // Fr point; // can be computed
    G1 eval;
    // std::vector<Fr> psi; //need to be deleted
    G1 w1;
    // G2 w2; //need to be deleted.

    KZGProof() = default;
    KZGProof& operator=(const KZGProof &other) = default;
    KZGProof(const KZGProof &other) = default;
    KZGProof(KZGProof &&other) = default;
    KZGProof(
        G1 &&eval,
        G1 &&w1):

    // point(std::move(point)),
    eval(std::move(eval)),
    w1(std::move(w1))
    {};

    size_t G1_size() const
    {
        return sizeof(w1) + sizeof(eval);
    }

    size_t size_in_bits() const
    {
       return (1 + w1.size_in_bits() + eval.size_in_bits());
    }

    void print_size() const
    {
        libff::print_indent(); printf("* G1 elements in Witness: %zu\n", this->G1_size());
        libff::print_indent(); printf("* Witness size in bits: %zu\n", this->size_in_bits());
    }

    bool operator==(const KZGProof &other) const;
    friend std::ostream& operator<< (std::ostream &out, const KZGProof &wit);
    friend std::istream& operator>> (std::istream &in, KZGProof &wit);
};

class CredProof {
public:
  // round 1
  G1 _A, _K;

  // round 2
  G1 _pi_tilde; // product_(j \in _I) pi_j
  GT _T_1, _T_2;
  

  // round 3
  Fr _mu;   // r_1 + r_2 * x
  Fr _r_x;  // r_3 * x + r_4 * x^2
  Fr _t_tilde;  // t(x)
  IPAProof _pi_ipa;   // 

  G1 _commit_fzy_1, _commit_fzy_2;
  G1 _commit_fzy, _commit_fg, _commit_fh;
  KZGProof _pi_kzg_fzy_1, _pi_kzg_fzy_2;
  KZGProof _pi_kzg_fzy;
  KZGProof _pi_kzg_fg;
  KZGProof _pi_kzg_fh;

  inline const size_t size_base() const {
    return 3 * sizeof(Fr) + 2 * sizeof(GT) + 3 * sizeof(G1) + _pi_ipa.size();
  }
  inline const size_t size_improved() const {
    return size_base() + 3 * sizeof(G1) + 3 * 2 * sizeof(G1); // 3 commitments, 3 times (2 G1 in a kzg proof)
  }
};

// 使用位运算判断是否只有一个位为1
bool isPowerOfTwo(size_t n);

Fr generate_random_sha256(const std::string& pre_image);

// Fr generate_random(string s);

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

template<typename T>
T sum_up(const std::vector<T>& v) {
  T sum = T::zero();
  for(auto item: v) {
    sum = sum + item;
  }
  return sum;
}

std::vector<Fr> multi_exponentiation_nlogn(const std::vector<Fr>& randoms, const std::vector<Fr>& randoms_inv, const size_t recursion_time);

std::vector<Fr> multi_exponentiation_n(const std::vector<Fr>& randoms, const size_t recursion_time);

/*
  poly_1: a_j x^j + ... + a_l x^l, base_1 = l, 默认小端序
  poly_2: a_k x^k + ... + a_m x^m, base_2 = m, 默认小端序

  result: a_n x^n + ... + a_0 x^0, base = 0, len = n, 默认大端序

  大端模式：高位字节存放在低地址中，低位字节存放在高地址中。最直观的字节序。
  小端模式：高位字节存放在高地址中，低位字节存放在低地址中。

  默认 poly_1 和 poly_2 的项都是连续的
*/
std::vector<Fr> poly_merge(
  const std::vector<Fr>& poly_1, const size_t base_1, \
  const std::vector<Fr>& poly_2, const size_t base_2, \
  const size_t n \
);

void test_poly_merge();

bool vector_all_zero(const std::vector<Fr>& v);

class AgendaItem {
  string label;
  long long enter_time;
  long long enter_cpu_time;
  long long leave_time;
  long long leave_cpu_time;

  size_t times_enter;
  size_t times_leave;
public:
  AgendaItem() {};
  AgendaItem(const string& label);
  void mark_enter();
  void mark_leave();
  void print(long long& start_time, long long& start_cpu_time);
  inline double duration() {
#ifdef CRED_DEBUG
    assert(times_enter == times_leave);
#endif
    long long time_from_last = (leave_time - enter_time) / this->times_enter;
    return time_from_last * 1e-9;
  }
};

class Agenda {
  string mode;  // base or improved
  size_t N;     // max bits
  size_t D;     // total bits
  size_t num;   // 
  size_t bit_len;

  size_t ipa_size;
  size_t proof_size;

  std::vector<string> label_records;
  std::map<string, AgendaItem> items;
  std::map<string, size_t> v_sizes;

  long long start_time;
  long long start_cpu_time;
public:
  Agenda() {};
  Agenda(const string mode, const size_t N, const size_t D, const size_t num, const size_t bit_len);
  void create_item(const string& label);
  void mark_item_end(const string& label);
  void mark_mem_usage(const string& label) {
    struct proc_t usage;
    look_up_our_self(&usage);
    if(v_sizes.find(label) == v_sizes.end()) {
      v_sizes[label] = usage.vsize;
    } else {
      v_sizes[label] = v_sizes[label] > usage.vsize ? v_sizes[label] : usage.vsize;
    }
  }
  void mark_proof_size(const CredProof& pi) {
    if(mode == "base") {
      this->proof_size = pi.size_base();
    } else {
      this->proof_size = pi.size_improved();
    }
  }
  void mark_ipa_size(const CredProof& pi) {
    assert(mode == "base");
    this->ipa_size = pi._pi_ipa.size();
  }

  string to_string_rough();
  string ipa_record();
  string bullet_record();
  void print();
  void write_file(const char* filename);
};

bool write_file(const string& file_path, const string& content);


// Agenda agenda;

};
#endif