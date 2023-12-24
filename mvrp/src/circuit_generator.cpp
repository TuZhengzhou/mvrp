#include <iostream>
#include <vector>
#include "libff/common/utils.hpp"
#include "libff/algebra/curves/public_params.hpp"
#include "basic_types.hpp"
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

CircuitParaGenerator::CircuitParaGenerator(const size_t num_multi_gates, const size_t num_commit_input) {

  this->num_multi_gates_ = num_multi_gates;
  this->num_commit_input_ = num_commit_input;
  this->num_constraints_ = num_multi_gates + num_commit_input - 1;

  size_t i;
  a_L.resize(this->num_multi_gates_);
  a_R.resize(this->num_multi_gates_);
  a_O.resize(this->num_multi_gates_);

  /** set random initial values*/
  a_L[0] = Fr::random_element();
  for(i = 0; i < this->num_multi_gates_; i++) {
    a_R[i] = Fr::random_element();
  }

  /** evaluate the circuit */
  for(i = 0; i < this->num_multi_gates_; i++) {
    a_O[i] = a_L[i] * a_R[i];
    if(i < this->num_multi_gates_-1)
      a_L[i+1] = a_O[i];
  }

  /** set the value of v */
  v.resize(this->num_commit_input_);
  for(i = 0; i < this->num_commit_input_; i++) {
    v[i] = a_R[i];
  }
  
  /** generate the constraint matrix */
  constraints = vector<LinearCombination>(this->num_constraints_, LinearCombination());

  for(i = 0; i < this->num_commit_input_; i++) {
    constraints[i] += LinearCombination(Variable(VType::MultiplierRight, i), Fr(1));
    constraints[i] += LinearCombination(Variable(VType::Committed, i), Fr(1));
  }
  for(i = 0; i < this->num_multi_gates_ - 1; i++) {
    constraints[this->num_commit_input_ + i] += LinearCombination(Variable(VType::MultiplierLeft, i+1), Fr( 1));
    constraints[this->num_commit_input_ + i] += LinearCombination(Variable(VType::MultiplierOutput, i), Fr(-1));
  }

#ifdef CRED_DEBUG

  assert(check());
#endif
}

bool CircuitParaGenerator::check() {
  bool ret = true;

  Fr z = Fr::random_element();

  for (auto lc: this->constraints) {
    Fr sum = Fr(0);
    cout << lc.toString() << endl;
    for (auto term: lc.terms) {
      switch (term.first.type)
      {
      case VType::MultiplierLeft:
        sum += a_L[term.first.index] * term.second;
        break;
      
      case VType::MultiplierRight:
        sum += a_R[term.first.index] * term.second;
        break;
      
      case VType::MultiplierOutput:
        sum += a_O[term.first.index] * term.second;
        break;

      case VType::Committed:
        sum -= v[term.first.index] * term.second;
        break;
      
      case VType::One:
        sum -= term.second;
        break;
      
      default:
        break;
      }
    }
    ret &= (sum == Fr(0));
  }
  return ret;
}



}; // namespace cred