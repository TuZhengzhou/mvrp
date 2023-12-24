#include <vector>
#include <ostream>
#include "ipa.hpp"
#include "basic_types.hpp"
#include "math_operations.hpp"

namespace cred {
Fr IPAProveSystem::generate_random(const Fr& pre_random, const G1& L, const G1& R) {
  std::string str;
  libff::bit_vector bits;

  // str  = pre_random.as_bigint().data.toString();
  str += L.coord[2].toString(10);
  str += R.coord[2].toString(10);

  // return cred::generate_random(str);
  return cred::generate_random_sha256(str);
}

G1 IPAProveSystem::ipa_hash(
  const std::vector<G1>& g_vec, const std::vector<G1>& h_vec, \
  const std::vector<Fr>& al_vec, const std::vector<Fr>& ar_vec, \
  const std::vector<Fr>& bl_vec, const std::vector<Fr>& br_vec, \
  const Fr& c) 
{
  size_t n      = g_vec.size();
  size_t n_apos = al_vec.size();
#ifdef CRED_DEBUG
  assert(h_vec.size() == n);
  assert(ar_vec.size() == n_apos);
  assert(bl_vec.size() == n_apos);
  assert(br_vec.size() == n_apos);
  assert(n_apos * (size_t)2 == n);
#endif 
  
  G1 ret        = c * this->_u;
  for(size_t i = 0; i < n_apos; i++) {
    ret = ret + al_vec[i] * g_vec[i] + ar_vec[i] * g_vec[n_apos+i];
    ret = ret + bl_vec[i] * h_vec[i] + br_vec[i] * h_vec[n_apos+i];
  }
  return ret;
}

G1 IPAProveSystem::ipa_hash(
  const std::vector<G1>& g_vec, const std::vector<G1>& h_vec, \
  const Fr& al, const Fr& ar, const Fr& bl, const Fr& br, const Fr& c
  ) 
{
  return ipa_hash(g_vec, h_vec, std::vector<Fr>(1, al), std::vector<Fr>(1, ar), std::vector<Fr>(1, bl), std::vector<Fr>(1, br), c);
}

G1 IPAProveSystem::ipa_hash(
  std::vector<G1>::const_iterator g_iter , std::vector<G1>::const_iterator h_iter, \
  std::vector<Fr>::const_iterator al_iter, std::vector<Fr>::const_iterator ar_iter, \
  std::vector<Fr>::const_iterator bl_iter, std::vector<Fr>::const_iterator br_iter, \
  const Fr& c, const size_t n_apos) 
{
  G1 ret        = c * this->_u;
  for(size_t i = 0; i < n_apos; i++) {
    ret = ret + (*(al_iter+i)) * (*(g_iter+       i));
    ret = ret + (*(ar_iter+i)) * (*(g_iter+n_apos+i));
    ret = ret + (*(bl_iter+i)) * (*(h_iter+       i));
    ret = ret + (*(br_iter+i)) * (*(h_iter+n_apos+i));
  }
  return ret;
}

G1 IPAProveSystem::ipa_hash(
    std::vector<G1>::const_iterator g_iter, std::vector<G1>::const_iterator h_iter, \
    std::vector<Fr>::const_iterator a_iter, std::vector<Fr>::const_iterator b_iter, \
    const Fr& c, const size_t n_apos, \
    char mode_a, char mode_b
  )
{
#ifdef CRED_DEBUG
  assert((mode_a == 'l' && mode_b == 'r') || (mode_a == 'r' && mode_b == 'l'));
#endif
  G1 ret        = c * this->_u;
  if(mode_a == 'l'){
    for(size_t i = 0; i < n_apos; i++) {
      ret = ret + (*(a_iter+i)) * (*(g_iter+       i));
      ret = ret + (*(b_iter+i)) * (*(h_iter+n_apos+i));
    }
  } else {
    for(size_t i = 0; i < n_apos; i++) {
      ret = ret + (*(a_iter+i)) * (*(g_iter+n_apos+i));
      ret = ret + (*(b_iter+i)) * (*(h_iter+       i));
    }
  }

  return ret;
}

IPAProveSystem::IPAProveSystem(G1& u){
   _u = u;
}

IPAProveSystem::IPAProveSystem(G1& P, G1& u, Fr& c){
   _P = P;
   _u = u;
   _c = c;
}

CIpaProofOneRecursion IPAProveSystem::IpaProveOneRecursion(
  const std::vector<G1>& g_vec, const std::vector<G1>& h_vec, \
  const std::vector<Fr>& a_vec, const std::vector<Fr>& b_vec, \
  const Fr& c
) {
  size_t n, n_apos;
  n = g_vec.size();
  n_apos = n / (size_t)2;
#ifdef CRED_DEBUG
  assert(h_vec.size() == n);
  assert(a_vec.size() == n);
  assert(b_vec.size() == n);
  assert(n_apos * (size_t)2 == n);
#endif

  CIpaProofOneRecursion pi;

  Fr cl, cr;
  std::vector<Fr> zero_vec = std::vector<Fr>(n_apos, Fr::zero());
  cl     = inner_product(a_vec.begin(), b_vec.begin() + n_apos, n_apos);
  cr     = inner_product(a_vec.begin() + n_apos, b_vec.begin(), n_apos);

  G1 L = ipa_hash(g_vec.begin()         , h_vec.begin(), \
                  zero_vec.begin()      , a_vec.begin(),\
                  b_vec.begin() + n_apos, zero_vec.begin(),\
                  cl, n_apos);
  G1 R = ipa_hash(g_vec.begin()         , h_vec.begin(), \
                  a_vec.begin() + n_apos, zero_vec.begin(), \
                  zero_vec.begin()      , b_vec.begin(), \
                  cr, n_apos);

  Fr x;
  x = generate_random(Fr::zero(), L, R);

  std::vector<Fr> a_apos, b_apos;
  a_apos = vector_add(numerical_mult(x, a_vec.begin(), n_apos), numerical_mult(x.inverse(), a_vec.begin()+n_apos, n_apos));
  b_apos = vector_add(numerical_mult(x.inverse(), b_vec.begin(), n_apos), numerical_mult(x, b_vec.begin()+n_apos, n_apos));

  pi.L = L;
  pi.R = R;

  pi.a_vec = a_apos;
  pi.b_vec = b_apos;

  return pi;
}


bool IPAProveSystem::IpaVerifyOneRecursion(
  const CIpaProofOneRecursion& pi, const std::vector<G1>& g_vec, const std::vector<G1>& h_vec
) {
  Fr c = inner_product(pi.a_vec, pi.b_vec);
  Fr x = generate_random(Fr::zero(), pi.L, pi.R);
  G1 P_apos = (x ^ 2) * pi.L + (x.inverse() ^ 2) * pi.R + (this->_P + this->_c * this->_u);

  G1 right  = ipa_hash(g_vec, h_vec, \
                      numerical_mult(x.inverse(), pi.a_vec), numerical_mult(x, pi.a_vec), \
                      numerical_mult(x, pi.b_vec), numerical_mult(x.inverse(), pi.b_vec), \
                      c);
  assert(P_apos == right);
  return P_apos == right;
} 

bool IPAOneRecursionTest() {
  size_t n  = (size_t)8;
  
  G1 u = Fr::random_element() * G1::one();
  std::vector<G1> g_vec(n, Fr::random_element() * G1::one());
  std::vector<G1> h_vec(n, Fr::random_element() * G1::one());
  std::vector<Fr> a_vec(n, Fr::random_element());
  std::vector<Fr> b_vec(n, Fr::random_element());
  Fr c = inner_product(a_vec, b_vec);

  G1 P = G1::zero();
  for(size_t i = 0; i < n; i++) {
    P = P + a_vec[i] * g_vec[i] + b_vec[i] * h_vec[i];
  }

  IPAProveSystem ipa_sys(P, u, c);
  CIpaProofOneRecursion pi = ipa_sys.IpaProveOneRecursion(g_vec, h_vec, a_vec, b_vec, c);
  bool result = ipa_sys.IpaVerifyOneRecursion(pi, g_vec, h_vec);
  assert(result == true);
  cout << "IPAOneRecursionTest result = " << result << endl;

  return result;
}

CIpaProof IPAProveSystem::IpaProve(
  const std::vector<G1>& g_vec, const std::vector<G1>& h_vec, \
  const std::vector<Fr>& a_vec, const std::vector<Fr>& b_vec, \
  const Fr& c) 
{
  size_t n = g_vec.size();
#ifdef CRED_DEBUG
  assert(h_vec.size() == n);
  assert(a_vec.size() == n);
  assert(b_vec.size() == n);

  assert(is_power_of_two(n));
#endif

  CIpaProof pi;
  std::vector<G1> g_apos = g_vec;
  std::vector<G1> h_apos = h_vec;
  std::vector<Fr> a_apos = a_vec;
  std::vector<Fr> b_apos = b_vec;

  size_t size, i, j, recursion_time;
  Fr x, pre_x, x_inv;
  std::vector<G1> g_apos_tmp, h_apos_tmp;
  std::vector<Fr> a_apos_tmp, b_apos_tmp;
  CIpaProofOneRecursion pi_one_recursion;

  // recursion_time = log2(n);
  pre_x = Fr::zero();
  size  = n / 2;

  while(a_apos.size() > (size_t)1) {

    pi_one_recursion = IpaProveOneRecursion(g_apos, h_apos, a_apos, b_apos, c);
    pi.L_vec.push_back(pi_one_recursion.L);
    pi.R_vec.push_back(pi_one_recursion.R);

    x    = generate_random(pre_x, pi_one_recursion.L, pi_one_recursion.R);
    x_inv= x.inverse();

    g_apos_tmp = std::vector<G1>(size);
    h_apos_tmp = std::vector<G1>(size);
    for(i = 0; i < size; i++) {
      g_apos_tmp[i] = (x_inv * g_apos[i]) + (x * g_apos[size+i]);
      h_apos_tmp[i] = (x * h_apos[i]) + (x_inv * h_apos[size+i]);
    }
    g_apos = g_apos_tmp;
    h_apos = h_apos_tmp;

    a_apos_tmp = std::vector<Fr>(size);
    b_apos_tmp = std::vector<Fr>(size);
    for(i = 0; i < size; i++) {
      a_apos_tmp[i] = (x * a_apos[i]) + (x_inv * a_apos[size+i]);
      b_apos_tmp[i] = (x_inv * b_apos[i]) + (x * b_apos[size+i]);
    }
    a_apos = a_apos_tmp;
    b_apos = b_apos_tmp;

    pre_x = x;
    size = size / 2;
  }
  pi.a = pi_one_recursion.a_vec[0];
  pi.b = pi_one_recursion.b_vec[0];
  return pi;
}

bool IPAProveSystem::IpaVerify(const CIpaProof& pi, const std::vector<G1>& g_vec, const std::vector<G1>& h_vec) {

  size_t i;
  size_t recursion_time = pi.L_vec.size();

  std::vector<Fr> randoms = std::vector<Fr>(recursion_time, Fr::zero());
  std::vector<Fr> randoms_inv = std::vector<Fr>(recursion_time, Fr::zero());
  libff::enter_block("IpaVerify::generate_random");
  randoms[0] = generate_random(Fr::zero(), pi.L_vec[0], pi.R_vec[0]);
  for(i = 1; i < recursion_time; i++) {
    randoms[i] = generate_random(randoms[i-1], pi.L_vec[i], pi.R_vec[i]);
  }
  libff::leave_block("IpaVerify::generate_random");
  for(i = 0; i < recursion_time; i++) {
    randoms_inv[i] = randoms[i].inverse();
  }

  G1 P_apos = (this->_P + this->_c * this->_u);
  for(size_t i = 0; i < recursion_time; i++) {
    P_apos = P_apos + (randoms[i] * randoms[i]) * pi.L_vec[i] + (randoms_inv[i] * randoms[i]) * pi.R_vec[i];
  }

  Fr c = pi.a * pi.b;

  std::vector<Fr> ss, ss_inv;
  ss = multi_exponentiation_n(randoms, recursion_time);
  ss_inv.resize(recursion_time);
  for(i = 0; i < recursion_time; i++) {
    ss_inv[i] = ss[i].inverse();
  }

  G1 g, h, right;
  g = G1::zero();
  h = G1::zero();
  for(i = 0; i < recursion_time; i++) {
    g = g + ss[i] * g_vec[i];
    h = h + ss[i].inverse() * h_vec[i];
  }

  right = pi.a * g + pi.b * h + c * this->_u;

  return P_apos == right;
}

bool IPAProveSystem::IpaMultiExpVerify(const CIpaProof& pi, const std::vector<G1>& g_vec, const std::vector<G1>& h_vec){
  size_t recursion_time = pi.L_vec.size();

  std::vector<Fr> randoms     = std::vector<Fr>(recursion_time, Fr::zero());
  std::vector<Fr> randoms_inv = std::vector<Fr>(recursion_time);
  libff::enter_block("IpaMultiExpVerify::generate_random");
  randoms[0] = generate_random(Fr::zero(), pi.L_vec[0], pi.R_vec[0]);
  for(size_t i = 1; i < recursion_time; i++) {
    randoms[i] = generate_random(randoms[i-1], pi.L_vec[i], pi.R_vec[i]);
  }
  for(size_t i = 0; i < recursion_time; i++) {
    randoms_inv[i] = randoms[i].inverse();
  }
  libff::leave_block("IpaMultiExpVerify::generate_random");

  size_t i;
  libff::enter_block("IpaMultiExpVerify::multi_exponentiation_n");
  std::vector<Fr> ss = multi_exponentiation_n(randoms, recursion_time);
  libff::leave_block("IpaMultiExpVerify::multi_exponentiation_n");
  
#ifdef CRED_DEBUG
  libff::enter_block("IpaMultiExpVerify::multi_exponentiation_nlogn");
  std::vector<Fr> ss1 = multi_exponentiation_nlogn(randoms, randoms_inv, recursion_time);
  libff::leave_block("IpaMultiExpVerify::multi_exponentiation_nlogn");

  std::vector<Fr> sub = vector_sub(ss, ss1);
  assert(vector_all_zero(sub));
#endif

  libff::enter_block("IpaMultiExpVerify::compute g");
  G1 g = G1::zero();
  G1 h = G1::zero();
  size_t len = g_vec.size();
  for(i = 0; i < len; i++) {
    g = g + ss[i] * g_vec[i];
  }
  libff::leave_block("IpaMultiExpVerify::compute g");

  libff::enter_block("IpaMultiExpVerify::compute h");
  for(i = 0; i < len; i++) {
    h = h + ss[i].inverse() * h_vec[i];
  }
  libff::leave_block("IpaMultiExpVerify::compute h");

  libff::enter_block("IpaMultiExpVerify::compute right, P_apos");
  G1 right = pi.a * g + pi.b * h + (pi.a * pi.b) * this->_u;

  G1 P_apos = (this->_P + this->_c * this->_u);
  for(size_t i = 0; i < recursion_time; i++) {
    P_apos = P_apos + (randoms[i] * randoms[i]) * pi.L_vec[i] + (randoms_inv[i] * randoms_inv[i]) * pi.R_vec[i];
  }
  libff::leave_block("IpaMultiExpVerify::compute right, P_apos");


  bool result = P_apos == right;
  return result;
}

bool IPAProveSystem::IpaMultiExpVerify(const CIpaProof& pi, const G1& g, const G1& h){
  
  size_t recursion_time = pi.L_vec.size();

  std::vector<Fr> randoms     = std::vector<Fr>(recursion_time, Fr::zero());
  std::vector<Fr> randoms_inv = std::vector<Fr>(recursion_time);

  libff::enter_block("IpaMultiExpVerify::generate_random");
  randoms[0] = generate_random(Fr::zero(), pi.L_vec[0], pi.R_vec[0]);
  for(size_t i = 1; i < recursion_time; i++) {
    randoms[i] = generate_random(randoms[i-1], pi.L_vec[i], pi.R_vec[i]);
  }
  for(size_t i = 0; i < recursion_time; i++) {
    randoms_inv[i] = randoms[i].inverse();
  }
  libff::leave_block("IpaMultiExpVerify::generate_random");


  G1 right = pi.a * g + pi.b * h + (pi.a * pi.b) * this->_u;

  G1 P_apos = (this->_P + this->_c * this->_u);
  for(size_t i = 0; i < recursion_time; i++) {
    P_apos = P_apos + (randoms[i] ^ 2) * pi.L_vec[i] + (randoms_inv[i] ^ 2) * pi.R_vec[i];
  }

  bool result = (P_apos == right);
  return result;
}

bool IPATest(size_t n) {
  
  G1 u = Fr::random_element() * G1::one();
  std::vector<G1> g_vec(n, Fr::random_element() * G1::one());
  std::vector<G1> h_vec(n, Fr::random_element() * G1::one());
  std::vector<Fr> a_vec(n, Fr::random_element());
  std::vector<Fr> b_vec(n, Fr::random_element());
  Fr c = inner_product(a_vec, b_vec);

  G1 P = G1::zero();
  for(size_t i = 0; i < n; i++) {
    P = P + a_vec[i] * g_vec[i] + b_vec[i] * h_vec[i];
  }

  bool result;

  IPAProveSystem ipa_sys(P, u, c);
  CIpaProof pi = ipa_sys.IpaProve(g_vec, h_vec, a_vec, b_vec, c);

  libff::enter_block("IpaVerify");
  result = ipa_sys.IpaVerify(pi, g_vec, h_vec);
  libff::leave_block("IpaVerify");

  libff::enter_block("IpaMultiExpVerify");
  result &= ipa_sys.IpaMultiExpVerify(pi, g_vec, h_vec);
  libff::leave_block("IpaMultiExpVerify");
  assert(result == true);

  cout << "IPATest result = " << result << endl;

  return result;
}
};  // cred