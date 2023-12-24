#include <openssl/sha.h>
#include "math_operations.hpp"

namespace cred {

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

std::vector<Fr> shift_and_reverse(const std::vector<Fr>& poly, const size_t shift) {
  size_t n = poly.size() + shift - 1;
  std::vector<Fr> ret(n + 1, Fr::zero());
  for(size_t i = 0; i < poly.size(); i++) {
    ret[n - (i + shift)] = poly[i];
  }
  return ret;
}
} // namespace cred