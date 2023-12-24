#include <iterator>
#include "structs.hpp"
#include "basic_types.hpp"

using namespace cred;

bool test_poly_merge() {
  size_t n1, n2, n;
  n1 = (size_t)3;
  n2 = (size_t)5;
  n  = (size_t)9;

  Fr zero = Fr::zero();
  Fr one  = zero + Fr::one();
  Fr two  = one  + Fr::one();
  Fr three  = two + Fr::one();
  Fr four  = three + Fr::one();
  Fr five  = four + Fr::one();
  Fr six  = five + Fr::one();
  Fr seven  = six + Fr::one();
  Fr eight  = seven + Fr::one();
  Fr nine  = eight + Fr::one();

  size_t base_1, base_2;
  base_1 = 1;
  base_2 = 5;
  std::vector<Fr> poly_1 = {one, two, three};
  std::vector<Fr> poly_2 = {five, six, seven, eight, nine};
  std::vector<Fr> poly = {nine, eight, seven, six, five, zero, three, two, one, zero};
  std::vector<Fr> poly_apos = poly_merge(poly_1, base_1, poly_2, base_2, n);

  std::ostream_iterator<Fr>   outFr(cout, ", ");

  libff::print_header("poly");
  copy(poly.begin(), poly.end(), outFr);
  cout << endl;

  libff::print_header("poly_apos");
  copy(poly_apos.begin(), poly_apos.end(), outFr);
  cout << endl;

  bool result = vector_all_zero(vector_sub(poly, poly_apos));

  return result;
}

int main() {
  libff::default_ec_pp::init_public_params();
  bool result = test_poly_merge();

  if (true == result)
    printf("test_poly_merge pass\n");
  return (int)(!result);
}