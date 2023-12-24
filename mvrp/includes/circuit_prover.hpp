#ifndef CRED_PROVER_CIRCUIT_HPP
#define CRED_PROVER_CIRCUIT_HPP
#include <vector>
#include <map>
#include <tuple>
#include "libff/algebra/curves/public_params.hpp"
#include "libff/algebra/fields/field_utils.hpp"
#include "libff/common/utils.hpp"
#include "structs.hpp"
#include "ipa.hpp"
#include "linear_combination.hpp"

namespace cred {

class ProverCircuit{
public:
  /** 公开输入 */
  size_t _N;    // 由 srs 推断

  size_t _num_constraints;    // constraints num, 由 WL 行数推断
  size_t _num_multi_gate;    // multi gates num, 由 WL 列数推断
  size_t _num_commit_input;   // commit input num, 由 WO 列数推断

  CCredSRS _srs; // srs, 外部输入
  G1 _V;    // commitment C, 初始化时计算得到
  vector<LinearCombination> _constraints; // 外部输入
  map<size_t, size_t> _map; // 外部输入

  /** 秘密输入 */
  Fr _r;    // random r for commitment C, 初始化时选取
  CSet _S;   // set S, 外部输入
  vector<Fr> _a_L; // 外部输入
  vector<Fr> _a_R; // 外部输入
  vector<Fr> _a_O; // 外部输入 n elements
  vector<Fr> _u;  // |I| elements (_num_commit_input elements) S + map + num_commit_input 推断

  /** 子证明系统和预计算 */
  IPAProveSystem _ipa_prove_sys; // 外部输入
  std::vector<G1> _point_proofs; // 初始化时计算得到
  G1 _beta_powers_sum_neg;

  ProverCircuit() {};
  ProverCircuit(const CCredSRS& srs, const vector<LinearCombination>& constraints, \
                const size_t num_multi_gate, const size_t num_commit_input, const map<size_t, size_t>& map, \
                const CSet& S, const vector<Fr>& aL, const vector<Fr>& aR, const vector<Fr>& aO, \
                const IPAProveSystem& ipa_sys);

  void flatten_constraints(const vector<LinearCombination>& constraints, const Fr& z, size_t num_multi_gate, vector<Fr>& wL, vector<Fr>& wR, vector<Fr>& wO);

  CCircuitProof prove(Agenda& agenda, const bool improved=true);
  CCircuitProof prove_base(Agenda& agenda);
  CCircuitProof prove_improved(Agenda& agenda);
  
private:
  inline const G1& point_proof(size_t idx);

  void commit_set();              // select random _r and generate commitment _C
  void point_proof_pre_compute(); // precompute pointproof for every set item
};

};
#endif