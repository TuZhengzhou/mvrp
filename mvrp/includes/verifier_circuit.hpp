#ifndef CRED_VERIFIER_CIRCUIT_HPP
#define CRED_VERIFIER_CIRCUIT_HPP
#include <vector>
#include <map>
#include <tuple>
#include "libff/algebra/curves/public_params.hpp"
#include "libff/algebra/fields/field_utils.hpp"
#include "structs.hpp"
#include "ipa.hpp"

namespace cred {

class VerifierCircuit{
public:
  /** 公开输入 */
  size_t _N;    // 由 srs 推断

  size_t _Q;    // constraints num, 由 WL 行数推断
  size_t _n;    // multi gates num, 由 WL 列数推断
  size_t _num_commit_input;   // commit input num, 由 WO 列数推断

  CredSRS _srs; // srs, 外部输入
  G1 _V;        // commitment C, 外部输入
  vector<vector<Fr>> _WL;   // 外部输入
  vector<vector<Fr>> _WR;   // 外部输入
  vector<vector<Fr>> _WO;   // 外部输入
  vector<vector<Fr>> _WU;   // 外部输入
  vector<Fr> _c;            // 外部输入
  vector<Fr> _I;            // 好像不需要这个东西, 可由 S + map + num_commit_input 推断
  map<size_t, size_t> _map; // 外部输入

  /** 子证明系统和预计算 */
  IPAProveSystem _ipa_prove_sys; // 外部输入
  G1 _beta_powers_sum_neg;

  VerifierCircuit() {};
  VerifierCircuit(const CredSRS& srs, const G1& V,\
                const vector<vector<Fr>>& WL, const vector<vector<Fr>>& WR, const vector<vector<Fr>>& WO, const vector<vector<Fr>>& WU, \
                vector<Fr> c, map<size_t, size_t> map, \
                const IPAProveSystem& ipa_sys);

  bool verify(const ProofCircuit& pi, Agenda& agenda, const bool improved = true);
  bool verify_base(const ProofCircuit& pi, Agenda& agenda);
  bool verify_improved(const ProofCircuit& pi, Agenda& agenda);

private:

};

  
};

#endif