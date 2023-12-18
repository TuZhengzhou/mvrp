#include <iostream>
#include <vector>
#include "libff/common/utils.hpp"
#include "libff/algebra/curves/public_params.hpp"
#include "structs.hpp"
#include "circuit_generator.hpp"

/**
 * this file supports generating circuits of the given pattern
 * l1 \
 *     * - o1(l2) \
 * r1 /            * - o2(l3) \
 *             r2 /            * - o3(l4)
 *                        r3 /            ...
 * 
 * let n be the number of multiplication gates.
 * let m be the number of committed inputs(v_i's).
 * 
 * to evaluate such a circuit, first we need to generate random inputs:
 *    l_1
 *    r_1, r_2,... r_n.
 * 
 * then we can evaluate the circuit, and output:
 *    [l_1, l_2,... l_n]
 *    [r_1, r_2,... r_n]
 *    [o_1, o_2,... o_n]
 * 
 * as for constrains of v_i's, we assume that r_i = v_i for i in [m]
 * so we have the first m constraints be like(for example, m = 3, n = 5):
 *     W_L              W_R              W_O             W_V          c
 *  0 0 0 0 0    |   1 0 0 0 0    |   0 0 0 0 0     |   1 0 0     |   0                       
 *  0 0 0 0 0    |   0 1 0 0 0    |   0 0 0 0 0     |   0 1 0     |   0                       
 *  0 0 0 0 0    |   0 0 1 0 0    |   0 0 0 0 0     |   0 0 1     |   0
 * 
 * we place the left n-1 constrains like:
 *     W_L              W_R              W_O             W_V          c
 *  0 1 0 0 0    |   0 0 0 0 0    |   -1 0 0 0 0     |   0 0 0     |   0                       
 *  0 0 1 0 0    |   0 0 0 0 0    |   0 -1 0 0 0     |   0 0 0     |   0                       
 *  0 0 0 1 0    |   0 0 0 0 0    |   0 0 -1 0 0     |   0 0 0     |   0                 
 *  0 0 0 0 1    |   0 0 0 0 0    |   0 0 0 -1 0     |   0 0 0     |   0                 
*/  

using namespace std;

namespace cred {

CircuitParaGenerator::CircuitParaGenerator(const size_t n, const size_t m) {
  // using cred::G1;
  // using cred::G2;
  // using cred::GT;
  // using cred::Fr;
  this->n = n;
  this->m = m;
  size_t i;
  a_L.resize(n);
  a_R.resize(n);
  a_O.resize(n);

  /** set random initial values*/
  a_L[0] = Fr::random_element();
  for(i = 0; i < n; i++) {
    a_R[i] = Fr::random_element();
  }

  /** evaluate the circuit */
  for(i = 0; i < n; i++) {
    a_O[i] = a_L[i] * a_R[i];
    if(i < n-1)
      a_L[i+1] = a_O[i];
  }

  /** set the value of v */
  v.resize(m);
  for(i = 0; i < m; i++) {
    v[i] = a_R[i];
  }
  
  /** generate the constraint matrix */
  W_L = vector<vector<Fr>>(n+m-1, vector<Fr>(n, Fr::zero()));
  W_R = vector<vector<Fr>>(n+m-1, vector<Fr>(n, Fr::zero()));
  W_O = vector<vector<Fr>>(n+m-1, vector<Fr>(n, Fr::zero()));
  W_V = vector<vector<Fr>>(n+m-1, vector<Fr>(m, Fr::zero()));
  c   = vector<Fr>(n+m-1, Fr::zero());
  for(i = 0; i < m; i++) {
    W_R[i][i] = Fr::one();
    W_V[i][i] = Fr::one();
  }
  for(i = 0; i < n-1; i++) {
    W_L[m+i][i+1] = Fr::one();
    W_O[m+i][i]   = - Fr::one();
  }
#ifdef CRED_DEBUG
  assert(check());

#endif
}

bool CircuitParaGenerator::check() {
  bool ret = true;
  ret &= vector_all_zero(vector_sub(a_O, hadmard_product(a_L, a_R)));
  for(size_t i = 0; i < n+m-1; i++) {
    Fr cal_result = inner_product(W_L[i], a_L) + inner_product(W_R[i], a_R) + inner_product(W_O[i], a_O) \
                    - inner_product(W_V[i], v) - c[i];
    ret &= (cal_result == Fr::zero());
  }
  return ret;
}



}; // namespace cred