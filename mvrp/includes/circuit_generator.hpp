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
#ifndef CRED_CIRCUIT_GENERATOR_HPP
#define CRED_CIRCUIT_GENERATOR_HPP

#include <vector>
#include "libff/algebra/curves/public_params.hpp"
#include "libff/common/default_types/ec_pp.hpp"
#include "structs.hpp"

using namespace std;

namespace cred {

class CircuitParaGenerator {
public:
  size_t n;
  size_t m;
  vector<Fr> a_L, a_R, a_O;
  vector<Fr> v;
  vector<vector<Fr>> W_L, W_R, W_O;
  vector<vector<Fr>> W_V;
  vector<Fr> c;
  CircuitParaGenerator(const size_t n, const size_t m);
  bool check();
};

}; // namespace cred

#endif; // CRED_CIRCUIT_GENERATOR_HPP