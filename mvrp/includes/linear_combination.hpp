#ifndef HPP_COMPILER_LINEAR_COMBINATION
#define HPP_COMPILER_LINEAR_COMBINATION

#include <vector>
#include <map>
#include <tuple>
#include <utility>
#include <iostream>
#include <string>
#include "../../mvrp/includes/basic_types.hpp"

using std::string;

// Assuming Scalar is a predefined type or class
// we can convert long to Fr by Fr(long)
using Scalar = cred::Fr;

// Function to convert Scalar to string
string to_string(const Scalar& s);

// Enum class for VType
enum class VType {
    Committed, 
    MultiplierLeft, 
    MultiplierRight, 
    MultiplierOutput, 
    One                 // For constant 1
};

std::string vTypeToString(VType type);


// Class for Variable
class Variable {
public:
    VType type;         // when type == VType::One, index is ignored
    size_t index;

    // Default constructor
    Variable() = default;

    // Constructor from VType and index
    Variable(VType type, size_t index) : type(type), index(index) {}

    // toString method
    string toString();

    // ... Other methods and operators as needed ...
};

// Class for LinearCombination
class LinearCombination {
public:
    std::vector<std::pair<Variable, Scalar>> terms;

    // Default constructor
    LinearCombination() = default;

    // Constructor from Variable
    LinearCombination(const Variable& var) {
        terms.emplace_back(var, Scalar(1));
    }

    // Constructor from Scalar
    LinearCombination(const Scalar& scalar) {
        terms.emplace_back(Variable(VType::One, 0), scalar);
    }

    // Constructor from a pair of Variable and Scalar
    LinearCombination(const Variable& var, const Scalar& scalar) {
        terms.emplace_back(var, scalar);
    }

    // Negation operator
    LinearCombination operator-() const {
        LinearCombination result;
        for (const auto& term : terms) {
            result.terms.emplace_back(term.first, -term.second);
        }
        return result;
    }

    // Add operation
    LinearCombination& operator+=(const LinearCombination& rhs) {
        terms.insert(terms.end(), rhs.terms.begin(), rhs.terms.end());
        return *this;
    }

    // Subtract operation
    LinearCombination& operator-=(const LinearCombination& rhs) {
        for (const auto& term : rhs.terms) {
            terms.emplace_back(term.first, -term.second);
        }
        return *this;
    }

    // Multiply operation
    LinearCombination& operator*=(const Scalar& rhs) {
        for (auto& term : terms) {
            term.second *= rhs;
        }
        return *this;
    }

    string toString();

    // ... Other methods and operators as needed ...
};

// Operator overloads for Variable
LinearCombination operator+(Variable var, const LinearCombination& lc);

LinearCombination operator-(Variable var, const LinearCombination& lc);

LinearCombination operator*(Variable var, const Scalar& scalar);

// Operator overloads for LinearCombination
LinearCombination operator+(const LinearCombination& lc1, const LinearCombination& lc2);

LinearCombination operator-(const LinearCombination& lc1, const LinearCombination& lc2);

// ... Other operator overloads as needed ...

std::ostream& operator<<(std::ostream& os, const Scalar& scalar);

std::ostream& operator<<(std::ostream& os, const Variable& var);

std::ostream& operator<<(std::ostream& os, const LinearCombination& lc);



#endif