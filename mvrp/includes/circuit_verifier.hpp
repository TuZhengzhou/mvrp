#ifndef CRED_VERIFIER_CIRCUIT_HPP
#define CRED_VERIFIER_CIRCUIT_HPP
#include <vector>
#include <map>
#include <tuple>
#include "libff/algebra/curves/public_params.hpp"
#include "libff/algebra/fields/field_utils.hpp"
#include "structs.hpp"
#include "ipa.hpp"
#include "linear_combination.hpp"

namespace cred {

class VerifierCircuit{
public:
  /** 公开输入 */
  size_t _N;    // 由 srs 推断

  size_t _num_constraints;    // constraints num
  size_t _num_multi_gate;     // number of multiplication gates
  size_t _num_commit_input;   // commit input num

  CCredSRS _srs; // srs, 外部输入
  vector<LinearCombination> _constraints; // 外部输入
  G1 _V;                    // commitment C, 外部输入
  vector<Fr> _I;            // 好像不需要这个东西, 可由 S + map + num_commit_input 推断
  map<size_t, size_t> _map; // 外部输入

  /** 子证明系统和预计算 */
  IPAProveSystem _ipa_prove_sys; // 外部输入
  G1 _beta_powers_sum_neg;

  VerifierCircuit() {};
  VerifierCircuit(const CCredSRS& srs, const vector<LinearCombination>& constraints, const G1& V,\
                const size_t num_multi_gate, const size_t num_commit_input, \
                const map<size_t, size_t>& map, const IPAProveSystem& ipa_sys);

  void flatten_constraints(const vector<LinearCombination>& constraints, const Fr& z, size_t num_multi_gate, \
                          Fr& wc, vector<Fr>& wL, vector<Fr>& wR, vector<Fr>& wO);

  bool verify(const CCircuitProof& pi, Agenda& agenda, const bool improved = true);
  bool verify_base(const CCircuitProof& pi, Agenda& agenda);
  bool verify_improved(const CCircuitProof& pi, Agenda& agenda);

private:

};

  
};

#endif