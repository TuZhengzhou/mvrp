#include <iostream>
#include <sstream>
#include <type_traits>
#include <libff/common/profiling.hpp>
#include <libff/algebra/fields/field_utils.hpp>
#include <libff/common/utils.hpp>
#include "basic_types.hpp"
#include "kzg.hpp"



namespace cred {

CKzgKey kzg_setup(const size_t t)
{
  libff::enter_block("Call to kzg_setup");

  /* Generate generator g, randomness a */
  libff::enter_block("Generator G, randomness A");

  G1 generator = G1::random_element();

  G2 generator2 = G2::random_element();  // why G2 together?: for the vefiryeval reduced_pairing computation

  Fr a = Fr::random_element();

  libff::leave_block("Generator G, randomness A");

  /* Generate t-SDH tuple : G1 */
  libff::enter_block("Generate t-SDH tuple: G1");

  Fr exp_a = Fr::one();

  std::vector<G1> g1tuple(t); //(t) 초기화 이후 emplace_back 연산 모두 배열 삽입 연산으로 바꿔줘야 한다.
  g1tuple[0] = generator; // t-SDH = (g1, ...)

  for(size_t i = 1; i < t; i++)
  {
    exp_a = exp_a * a;
    g1tuple[i] = exp_a * generator; //group element should be at right side always!! ALWAYS !!!!!
  }

  libff::leave_block("Generate t-SDH tuple: G1");

  /* Generate t-SDH tuple : G2 */
  libff::enter_block("Generate t-SDH tuple: G2");

  Fr exp_a2 = Fr::one();

  std::vector<G2> g2tuple(t);
  g2tuple[0] = generator2;

  for(size_t i = 1; i < t; i++)
  {
    exp_a2 = exp_a2 * a;
    g2tuple[i] = exp_a2 * generator2;
  }

  libff::leave_block("Generate t-SDH tuple: G2");

  /* Output as a commitment key */
  libff::leave_block("Call to kzg_setup");

  CKzgKey tuple = CKzgKey(std::move(g1tuple), std::move(g2tuple));
  return tuple;

}

G1 kzg_commit(const CKzgKey &ck, const std::vector<Fr> &poly, const size_t t)
{
  libff::enter_block("Call to kzg_commit");

  libff::enter_block("Commit at G1");

  G1 temp = G1::zero();
  G1 commit1 = G1::zero();

  for(size_t i = 1; i < t + 1; i++) {
    if(poly[t - i] == 0) {
      continue;
    }
    else {
      temp = poly[t - i] * ck.g1[i - 1];
      commit1 = temp + commit1;
    }
  }

  libff::leave_block("Commit at G1");

  libff::leave_block("Call to kzg_commit");

  return commit1;
}

Fr kzg_hash(const G1 &commit_a, const G1 &commit_b, const G1 &commit_c)
{
  libff::enter_block("Call to kzg_hash");
  //print <-> print_coordinates: differnet
  libff::enter_block("Extract Commit A, B, C's Coord[2]");
  std::string a_x = commit_a.coord[2].toString(10);
  std::string b_x = commit_b.coord[2].toString(10);
  std::string c_x = commit_c.coord[2].toString(10);
  a_x += (b_x + c_x);
  libff::leave_block("Extract Commit A, B, C's Coord[2]");

  return cred::generate_random_sha256(a_x);
}

Fr kzg_hash(const G1 &commit)
{
  libff::enter_block("Call to kzg_hash");
  //print <-> print_coordinates: differnet
  // libff::enter_block("Extract Commit A's Coord[2]");
  std::string a_x = commit.coord[2].toString(10);
  // libff::leave_block("Extract Commit A's Coord[2]");

  Fr random_point = cred::generate_random_sha256(a_x);
  libff::leave_block("Call to kzg_hash");
  return random_point;
}

CKzgProof kzg_prove(const CKzgKey &ck, std::vector<Fr> &poly, const Fr &point, size_t t)
{
  libff::enter_block("Call to kzg_prove");
  size_t i;

  /* Evaluate Polynomial */
  libff::enter_block("Evaluate Polynomial");

  Fr eval = kzg_evaluate(poly, point, t);

  libff::leave_block("Evaluate Polynomial");
  poly[t - 1] = poly[t - 1] - eval;
  std::vector<Fr> divisor(2);
  // divisor[0] = convert(1);
  divisor[0] = Fr::one();

  Fr minus = Fr::zero();
  minus = minus - point;
  divisor[1] = minus;

  // libff::leave_block("Compute Divisor[2]: stands for polynomial (x - point)");

  //division Algorithm.
  // libff::enter_block("Divide Algorithm: poly(x) - poly(i) / (x - i)");
  std::vector<Fr> psi(t);

  for(i = 0; i < t - 1; i++) {
    psi[i] = poly[i];
    poly[i] = poly[i] - (psi[i] * divisor[0]);
    poly[i + 1] = poly[i + 1] - psi[i] * divisor[1];
  }

  // libff::leave_block("Divide Algorithm: poly(x) - poly(i) / (x - i)");

  /* compute w = g ^ psi(a) */

  // libff::enter_block("Compute w = g ^ psi(a): G1");

  G1 temp1 = G1::zero();
  G1 w1 = G1::zero();

  for(size_t i = 2; i < t + 1; i++) {
    if(psi[t - i] == 0) {
      continue;
    }
    else {
      temp1 = psi[t - i] * ck.g1[i - 2];
      w1 = temp1 + w1;
    }
  }

  // libff::leave_block("Compute w = g ^ psi(a): G1");

  /* For Non-Interactive: Put evaluation as Group-1 Element */
  G1 eval_g1 = eval * ck.g1[0];
  
  libff::leave_block("Call to kzg_prove");

  /* Output as a CKzgProof */
  CKzgProof wit = CKzgProof(std::move(eval_g1), std::move(w1));
  wit.psi.resize(t - 1);
  for(i = 0; i < t - 1; i++) {
    wit.psi[i] = psi[i] - poly[t - 1];
  }

  return wit;
}

/** 将设置密钥、承诺、生成证明打包在一起 
 * exp_start: 多项式的最低项幂次
 * exp_span: 多项式最高项幂次-最低项幂次+1
 * poly: 多项式各项系数, 应当有 poly.size() = exp_span
 * 例如: 
 * 对于多项式 f(x) = 2 x^2 + 6 x^8
 * exp_start 为 2
 * exp_span 为 7(8-2+1)
//  * poly 为 [2,0,0,0,0,0,6]
*/
CKzgProof kzg_prove(const CCredSRS& srs, std::vector<Fr>& poly, size_t t, G1& commit) {
#ifdef CRED_DEBUG
  assert(poly.size() == t);
#endif
  std::vector<G1> g1tuple;
  std::vector<G2> g2tuple;
  g1tuple.insert(g1tuple.end(), srs.g1s.begin(), srs.g1s.begin() + t);
  g2tuple.insert(g2tuple.end(), srs.g2s.begin(), srs.g2s.begin() + (size_t)1);
  
  CKzgKey kzg_ck(g1tuple, g2tuple);

  Fr random;
  commit = kzg_commit(kzg_ck, poly, t);
  random = kzg_hash(commit);
  
  CKzgProof pi;
  pi = kzg_prove(kzg_ck, poly, random, t);

  GT beta_N_plus_2_GT = ReducedPairing(srs.get_g1_beta_exp(srs.N), srs.get_g2_beta_exp(2));

  if (pi.psi.size() > srs.N+1) {
    Fr coef_exp_N_plus_one;
    coef_exp_N_plus_one = pi.psi[pi.psi.size() - (srs.N+2)];
    pi.wt = ((beta_N_plus_2_GT) ^ (coef_exp_N_plus_one)) * (((srs.gt) ^ (coef_exp_N_plus_one * random)).unitary_inverse());
  } else {
    pi.wt = (srs.gt) ^ (Fr::zero());
  }
  pi.psi.clear();

  return pi;
}



Fr kzg_evaluate(const std::vector<Fr> &poly, const Fr &point, const size_t t)
{

  Fr eval = Fr::zero();
  Fr temp = Fr::one();

  for(size_t i = 1; i < t + 1; i++) {
  eval += poly[t - i] * temp;
  temp *= point;
  }

  return eval;
}

bool kzg_vfyeval(const CKzgKey& ck, const G1& commit, const CKzgProof &CKzgProof)
{
  libff::enter_block("Call to kzg_vfyeval");
  GT left1 = ReducedPairing(commit, ck.g2[0]); //either side does not matter.
  Fr zero = Fr::zero();

  G2 num = kzg_hash(commit) * ck.g2[0];

  GT right1 = ReducedPairing(CKzgProof.w1, ck.g2[1] - num);
  GT right2 = ReducedPairing(CKzgProof.eval, ck.g2[0]); //eval which side? doesnt matter.

  GT right = right1 * right2 * CKzgProof.wt;
  bool verifyresult;
  if (left1 == right) {
    verifyresult = true;
  } else {
    verifyresult = false;
  }
  libff::leave_block("Call to kzg_vfyeval");
  return verifyresult;
}

bool kzg_test() {
  libff::enter_block("Call to run_kzg");

  size_t t = 10000;

  /* Generate Polynomial to Commit: we need to put Convolution Poly. in this section */
  // 一个 t-1 阶多项式，有 t 个系数
  std::vector<Fr> poly(t);

  for(size_t i = 0; i < t; i++) {
    Fr random = Fr::random_element();
    poly[i] = random;
    // poly[i].print();
  }

  /* Generate t-SDH tuple, and select secret randomness t */

  libff::print_header("Generate Key: t-SDH Tuple");
  CKzgKey ck = kzg_setup(t);
  printf("\n"); libff::print_indent(); libff::print_mem("after setup");

  /* Commit Polynomial into Product: G1-element */

  libff::print_header("Commit Polynomial");
  G1 commit = kzg_commit(ck, poly, t);
  printf("\n"); libff::print_indent(); libff::print_mem("after commit");

  /* Generate Random Point for Evaluation */

  Fr point = kzg_hash(commit);
  // point.print();

  /* Generate Witness of the evaluation + Evaluate the Polynomial */

  libff::print_header("Create Witness");
  CKzgProof wit = kzg_prove(ck, poly, point, t);
  printf("\n"); libff::print_indent(); libff::print_mem("after create-CKzgProof");

  /* Verify evaluation */
  libff::print_header("Verify Evaluation of Polynomial");
  bool verifyresult = kzg_vfyeval(ck, commit, wit);

  if (verifyresult == true) {
    libff::print_header("VERIFICATION ACCEPT!!");
  } else {
    libff::print_header("VERIFICATION REJECT");
  }
  
  printf("\n"); libff::print_indent(); libff::print_mem("after vfyeval");

  poly.clear();
  point.clear();

  libff::leave_block("Call to run_kzg");

  return verifyresult;
}

} // namespace cred
