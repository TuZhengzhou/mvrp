#ifndef MATH_OPERATIONS_HPP
#define MATH_OPERATIONS_HPP

#include <vector>
#include "libff/common/profiling.hpp"
#include "libff/algebra/curves/public_params.hpp"
#include "libff/common/default_types/ec_pp.hpp"
#include "basic_types.hpp"

using std::vector;

namespace cred {

Fr inner_product(const std::vector<Fr> &v1, const std::vector<Fr> &v2);

// need to check if the relation of v1, v2, len
Fr inner_product(std::vector<Fr>::const_iterator v1_iter, std::vector<Fr>::const_iterator v2_iter, size_t len);

std::vector<Fr> hadmard_product(const std::vector<Fr> &v1, const std::vector<Fr> &v2);

std::vector<Fr> numerical_mult(const Fr &factor, const std::vector<Fr> &v);

std::vector<Fr> numerical_mult(const Fr &factor, std::vector<Fr>::const_iterator v_begin, size_t size);

std::vector<Fr> vector_add(const std::vector<Fr> &v1, const std::vector<Fr> &v2);

std::vector<Fr> vector_sub(const std::vector<Fr> &v1, const std::vector<Fr> &v2);

std::vector<Fr> vector_neg(const std::vector<Fr> &v);

/*
  input: x, max_exp
  output: 
    1) x^0, x^1, x^2, ..., x^max_exp if zero_exp = true
    2)      x^1, x^2, ..., x^max_exp if zero_exp = false

  zero_exp is set to true by default
*/
std::vector<Fr> vector_powers(const Fr &x, size_t max_exp, bool zero_exp = true);

std::vector<Fr> vec_matrix_mult(const vector<Fr>& vec, const vector<vector<Fr>>& matrx, const bool vec_left);

// a special case of inner_product where the v2 vector is a vector of 1
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

std::vector<Fr> shift_and_reverse(const std::vector<Fr>& poly, const size_t shift);



// void test_poly_merge();

bool vector_all_zero(const std::vector<Fr>& v);

} // namespace cred

#endif // MATH_OPERATIONS_HPP