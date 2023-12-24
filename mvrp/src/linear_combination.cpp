// 
#include "linear_combination.hpp"

// Function to convert Scalar to string
string to_string(const Scalar& s) {
    // return std::to_string(s);
    return std::to_string(s.as_ulong());
}

std::string vTypeToString(VType type) {
    switch (type) {
        case VType::Committed:        return "Committed";
        case VType::MultiplierLeft:   return "MultiplierLeft";
        case VType::MultiplierRight:  return "MultiplierRight";
        case VType::MultiplierOutput: return "MultiplierOutput";
        case VType::One:              return "One";
        default:                      return "Unknown";
    }
}

// Class for Variable
string Variable::toString() {
    if (type == VType::One) {
        return "1";
    } else {
        return vTypeToString(type) + "(" + std::to_string(index) + ")";
    }
}

string LinearCombination::toString() {
    string result = "";
    for (size_t i = 0; i < terms.size(); ++i) {
        if (i > 0) result += " + ";
        result += to_string(terms[i].second) + "*" + terms[i].first.toString();
    }
    return result;
}

// Operator overloads for Variable
LinearCombination operator+(Variable var, const LinearCombination& lc) {
    return LinearCombination(var, Scalar(1)) += lc;
}

LinearCombination operator-(Variable var, const LinearCombination& lc) {
    return LinearCombination(var, Scalar(1)) -= lc;
}

LinearCombination operator*(Variable var, const Scalar& scalar) {
    return LinearCombination(var, scalar);
}

// Operator overloads for LinearCombination
LinearCombination operator+(const LinearCombination& lc1, const LinearCombination& lc2) {
    return LinearCombination(lc1) += lc2;
}

LinearCombination operator-(const LinearCombination& lc1, const LinearCombination& lc2) {
    return LinearCombination(lc1) -= lc2;
}

// ... Other operator overloads as needed ...

std::ostream& operator<<(std::ostream& os, const Scalar& scalar) {
    os << to_string(scalar);
    return os;
}

std::ostream& operator<<(std::ostream& os, const Variable& var) {
    os << Variable(var).toString();
    return os;
}

std::ostream& operator<<(std::ostream& os, const LinearCombination& lc) {
    os << LinearCombination(lc).toString();
    return os;
}

