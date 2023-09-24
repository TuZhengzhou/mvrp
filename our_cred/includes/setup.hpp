#ifndef CRED_SETUP_HPP
#define CRED_SETUP_HPP
// #include "../../depends/libff/libff/algebra/curves/"
#include "libff/algebra/curves/bn128/bn128_pp.hpp"
#include "./structs.hpp"

namespace cred {
/*
choose:
  1) g1, g2 : generator of G1, G2
  2) β : random

calculate:
  1) 2N - 1 in G1
  2) N in G2
  3) 1 in GT
*/
CredSRS trustSetup(uint64_t N);
};
#endif