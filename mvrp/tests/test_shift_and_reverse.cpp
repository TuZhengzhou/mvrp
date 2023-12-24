#include "structs.hpp"

namespace cred {

void testReverseAndShift() {
    libff::default_ec_pp::init_public_params();
    
    // Test with a simple case
    std::vector<Fr> poly;
    poly.push_back(Fr(1));
    poly.push_back(Fr(2));
    poly.push_back(Fr(3));
    size_t base = 2;
    toString(poly);
    
    auto result = shift_and_reverse(poly, base);
    toString(result);

    // Expected result: {3, 2, 1, 0, 0}
    assert(result.size() == 5);
    assert(result[0] == Fr(3));
    assert(result[1] == Fr(2));
    assert(result[2] == Fr(1));
    assert(result[3] == Fr(0));
    assert(result[4] == Fr(0));

    std::cout << "Test passed.\n";
}

}

int main() {
    cred::testReverseAndShift();
    return 0;
}