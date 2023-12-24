#include "linear_combination.hpp"
#include <iostream>
#include <sstream>
#include <cassert>

int testVariable() {
    Variable v1(VType::MultiplierLeft, 2);
    assert(v1.toString() == "MultiplierLeft(2)");
    if (v1.toString() != "MultiplierLeft(2)") return 1;

    Variable v2(VType::One, 0);
    assert(v2.toString() == "1");
    if (v2.toString() != "1") return 1;

    std::cout << "Variable tests passed.\n";
    return 0;
}

int testLinearCombination() {
    Variable v1(VType::MultiplierLeft, 2);
    Scalar s1 = Scalar(3);

    LinearCombination lc1(v1, s1);

    assert(lc1.toString() == "3*MultiplierLeft(2)");
    
    if (lc1.toString() != "3*MultiplierLeft(2)") return 1;

    Variable v2(VType::One, 0);
    Scalar s2 = Scalar(5);
    LinearCombination lc2(s2);
    assert(lc2.toString() == "5*1");
    if (lc2.toString() != "5*1") return 1;

    std::cout << "LinearCombination tests passed.\n";
    return 0;
}

int testOperations() {
    Variable v1(VType::MultiplierLeft, 2);
    Scalar s1 = Scalar(3);
    LinearCombination lc1(v1, s1);

    Variable v2(VType::MultiplierRight, 3);
    Scalar s2 = Scalar(4);
    LinearCombination lc2(v2, s2);

    // Test addition
    LinearCombination lcAdd = lc1 + lc2;
    assert(lcAdd.toString() == "3*MultiplierLeft(2) + 4*MultiplierRight(3)");
    if (lcAdd.toString() != "3*MultiplierLeft(2) + 4*MultiplierRight(3)") return 1;

    // Test subtraction
    LinearCombination lcSub = lc1 - lc2;
    assert(lcSub.toString() == "3*MultiplierLeft(2) + 4891460686036598781*MultiplierRight(3)");
    if (lcSub.toString() != "3*MultiplierLeft(2) + 4891460686036598781*MultiplierRight(3)") return 1;

    // Test negation
    LinearCombination lcNeg = -lc1;
    assert(lcNeg.toString() == "4891460686036598782*MultiplierLeft(2)");
    if (lcNeg.toString() != "4891460686036598782*MultiplierLeft(2)") return 1;

    // Test multiplication
    lc1 *= 2;
    assert(lc1.toString() == "6*MultiplierLeft(2)");
    if (lc1.toString() != "6*MultiplierLeft(2)") return 1;

    std::cout << "Operations tests passed.\n";
    return 0;
}

int testOutput() {
    Variable v(VType::MultiplierLeft, 3);
    Scalar s = Scalar(5);
    LinearCombination lc(v, s);

    std::ostringstream os;
    os << lc;
    assert(os.str() == "5*MultiplierLeft(3)");
    if (os.str() != "5*MultiplierLeft(3)") return 1;

    std::cout << "Output tests passed.\n";
    return 0;
}

int main() {
    libff::default_ec_pp::init_public_params();

    if (testVariable() != 0) return 1;
    if (testLinearCombination() != 0) return 1;
    if (testOperations() != 0) return 1;
    if (testOutput() != 0) return 1;

    std::cout << "All tests passed!\n";
    return 0;
}
