#ifndef CRED_PROVER_HPP
#define CRED_PROVER_HPP
#include <vector>
#include <map>
#include <tuple>
#include "libff/algebra/curves/public_params.hpp"
#include "libff/algebra/fields/field_utils.hpp"
#include "libff/common/utils.hpp"
#include "structs.hpp"
#include "range.hpp"

namespace cred {

class Prover{
public:
  CCredSRS _srs; // srs
  CSet _S;   // set S
  G1 _C;    // commitment C
  Fr _r;    // random r for commitment C
  CRanges _ranges;
  IPAProveSystem _ipa_prove_sys;

  size_t _D;    // D = sum_(j in I) n_(j)
  size_t _N;    // _N = srs.N
  map<size_t, size_t> _part_sums; // sum_(i=1)^(j-1) n_i

  std::vector<G1> _point_proofs;

  Prover() {};
  Prover(const CCredSRS& srs, const CSet& S, const CRanges& ranges, const IPAProveSystem& ipa_sys);

  CRangeProof prove(Agenda& agenda, const bool improved=true);
  CRangeProof prove_base(Agenda& agenda);
  CRangeProof prove_improved(Agenda& agenda);
  
private:
  inline const G1& point_proof(size_t idx);

  void commit_set();              // select random _r and generate commitment _C
  void point_proof_pre_compute(); // precompute pointproof for every set item
  bool t0_check(const Fr& two, const vector<Fr>& one_vec, const vector<Fr>& z_powers_0, const vector<Fr>& y_powers_1);
  bool easy_check(const std::vector<Fr>& z_exps); // 直接给元素值, 然后检查
};

};
#endif