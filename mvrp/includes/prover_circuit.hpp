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

namespace cred {

class ProverCircuit{
public:
  /** 公开输入 */
  size_t _N;    // 由 srs 推断

  size_t _Q;    // constraints num, 由 WL 行数推断
  size_t _n;    // multi gates num, 由 WL 列数推断
  size_t _num_commit_input;   // commit input num, 由 WO 列数推断

  CredSRS _srs; // srs, 外部输入
  G1 _V;    // commitment C, 初始化时计算得到
  vector<vector<Fr>> _WL;   // 外部输入
  vector<vector<Fr>> _WR;   // 外部输入
  vector<vector<Fr>> _WO;   // 外部输入
  vector<vector<Fr>> _WU;   // 外部输入
  vector<Fr> _c;            // 外部输入
  vector<Fr> _I;            // 好像不需要这个东西, 可由 S + map + num_commit_input 推断
  map<size_t, size_t> _map; // 外部输入

  /** 秘密输入 */
  Fr _r;    // random r for commitment C, 初始化时选取
  SET _S;   // set S, 外部输入
  vector<Fr> _a_L; // 外部输入
  vector<Fr> _a_R; // 外部输入
  vector<Fr> _a_O; // 外部输入 n elements
  vector<Fr> _u;  // |I| elements (_num_commit_input elements) S + map + num_commit_input 推断

  /** 子证明系统和预计算 */
  IPAProveSystem _ipa_prove_sys; // 外部输入
  std::vector<G1> _point_proofs; // 初始化时计算得到
  G1 _beta_powers_sum_neg;

  ProverCircuit() {};
  ProverCircuit(const CredSRS& srs, \
                const vector<vector<Fr>>& WL, const vector<vector<Fr>>& WR, const vector<vector<Fr>>& WO, const vector<vector<Fr>>& WU, \
                vector<Fr> c, map<size_t, size_t> map, \
                const SET& S, vector<Fr> aL, vector<Fr> aR, vector<Fr> aO, \
                const IPAProveSystem& ipa_sys);

  ProofCircuit prove(Agenda& agenda, const bool improved=true);
  ProofCircuit prove_base(Agenda& agenda);
  ProofCircuit prove_improved(Agenda& agenda);
  
private:
  inline const G1& point_proof(size_t idx);

  void commit_set();              // select random _r and generate commitment _C
  void point_proof_pre_compute(); // precompute pointproof for every set item
};

};
#endif