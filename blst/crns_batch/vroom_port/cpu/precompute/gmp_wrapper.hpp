#ifndef GMP_WRAPPER_HPP
#define GMP_WRAPPER_HPP

#include <gmpxx.h>
#include <string>
#include <iostream>
#include <stdexcept>
#include <ctime>

class BigInt {
private:
    mpz_class value;

public:
    // Constructors
    BigInt() : value(0) {}
    BigInt(const mpz_class& val) : value(val) {}
    BigInt(const mpz_t val) : value(val) {}
    BigInt(long val) : value(val) {}
    BigInt(unsigned long val) : value(val) {}
    BigInt(int val) : value(val) {}
    BigInt(unsigned int val) : value(val) {}
    BigInt(const std::string& str, int base = 10) : value(str, base) {}
    BigInt(const char* str, int base = 10) : value(str, base) {}
    
    // Copy constructor
    BigInt(const BigInt& other) : value(other.value) {}
    
    // Assignment operator
    BigInt& operator=(const BigInt& other) {
        if (this != &other) {
            value = other.value;
        }
        return *this;
    }
    
    // Assignment from other types
    BigInt& operator=(const mpz_class& val) {
        value = val;
        return *this;
    }
    BigInt& operator=(long val) {
        value = val;
        return *this;
    }
    BigInt& operator=(const std::string& str) {
        value = mpz_class(str);
        return *this;
    }

    // Arithmetic operators
    BigInt operator+(const BigInt& other) const {
        return BigInt(value + other.value);
    }
    
    BigInt operator-(const BigInt& other) const {
        return BigInt(value - other.value);
    }
    
    BigInt operator*(const BigInt& other) const {
        return BigInt(value * other.value);
    }
    
    BigInt operator/(const BigInt& other) const {
        if (other.value == 0) {
            throw std::runtime_error("Division by zero");
        }
        return BigInt(value / other.value);
    }
    
    BigInt operator%(const BigInt& other) const {
        if (other.value == 0) {
            throw std::runtime_error("Modulo by zero");
        }
        return BigInt(value % other.value);
    }
    
    // Unary operators
    BigInt operator-() const {
        return BigInt(-value);
    }
    
    BigInt operator+() const {
        return *this;
    }
    
    // Compound assignment operators
    BigInt& operator+=(const BigInt& other) {
        value += other.value;
        return *this;
    }
    
    BigInt& operator-=(const BigInt& other) {
        value -= other.value;
        return *this;
    }
    
    BigInt& operator*=(const BigInt& other) {
        value *= other.value;
        return *this;
    }
    
    BigInt& operator/=(const BigInt& other) {
        if (other.value == 0) {
            throw std::runtime_error("Division by zero");
        }
        value /= other.value;
        return *this;
    }
    
    BigInt& operator%=(const BigInt& other) {
        if (other.value == 0) {
            throw std::runtime_error("Modulo by zero");
        }
        value %= other.value;
        return *this;
    }

    // Comparison operators
    bool operator<(const BigInt& other) const {
        return value < other.value;
    }
    
    bool operator>(const BigInt& other) const {
        return value > other.value;
    }
    
    bool operator<=(const BigInt& other) const {
        return value <= other.value;
    }
    
    bool operator>=(const BigInt& other) const {
        return value >= other.value;
    }
    
    bool operator==(const BigInt& other) const {
        return value == other.value;
    }
    
    bool operator!=(const BigInt& other) const {
        return value != other.value;
    }

    // Bitwise operators
    BigInt operator<<(int shift) const {
        return BigInt(value << shift);
    }
    
    BigInt operator>>(int shift) const {
        return BigInt(value >> shift);
    }
    
    BigInt operator&(const BigInt& other) const {
        return BigInt(value & other.value);
    }
    
    BigInt operator|(const BigInt& other) const {
        return BigInt(value | other.value);
    }
    
    BigInt operator^(const BigInt& other) const {
        return BigInt(value ^ other.value);
    }
    
    BigInt operator~() const {
        return BigInt(~value);
    }
    
    // Compound bitwise assignment operators
    BigInt& operator<<=(int shift) {
        value <<= shift;
        return *this;
    }
    
    BigInt& operator>>=(int shift) {
        value >>= shift;
        return *this;
    }
    
    BigInt& operator&=(const BigInt& other) {
        value &= other.value;
        return *this;
    }
    
    BigInt& operator|=(const BigInt& other) {
        value |= other.value;
        return *this;
    }
    
    BigInt& operator^=(const BigInt& other) {
        value ^= other.value;
        return *this;
    }

    // Power method with modular exponentiation
    // For negative exponents, computes modular inverse
    BigInt pow(const BigInt& exponent, const BigInt& modulus = BigInt(0)) const {
        if (modulus == 0) {
            // Regular exponentiation without modulus
            if (exponent < 0) {
                throw std::runtime_error("Negative exponent not supported for non-modular exponentiation");
            }
            mpz_class result;
            mpz_pow_ui(result.get_mpz_t(), value.get_mpz_t(), mpz_get_ui(exponent.value.get_mpz_t()));
            return BigInt(result);
        } else {
            // Modular exponentiation
            mpz_class result;
            if (exponent >= 0) {
                mpz_powm(result.get_mpz_t(), value.get_mpz_t(), 
                        exponent.value.get_mpz_t(), modulus.value.get_mpz_t());
            } else {
                // For negative exponents, compute modular inverse first
                BigInt abs_exp = -exponent;
                mpz_class inverse;
                if (mpz_invert(inverse.get_mpz_t(), value.get_mpz_t(), modulus.value.get_mpz_t()) == 0) {
                    throw std::runtime_error("Modular inverse does not exist");
                }
                mpz_powm(result.get_mpz_t(), inverse.get_mpz_t(), 
                        abs_exp.value.get_mpz_t(), modulus.value.get_mpz_t());
            }
            return BigInt(result);
        }
    }

    // Utility methods
    std::string to_string(int base = 10) const {
        return value.get_str(base);
    }
    
    long to_long() const {
        return value.get_si();
    }
    
    unsigned long to_ulong() const {
        return value.get_ui();
    }
    
    double to_double() const {
        return value.get_d();
    }
    
    int sign() const {
        return mpz_sgn(value.get_mpz_t());
    }
    
    size_t bit_length() const {
        return mpz_sizeinbase(value.get_mpz_t(), 2);
    }
    
    bool is_zero() const {
        return value == 0;
    }
    
    bool is_negative() const {
        return value < 0;
    }
    
    bool is_positive() const {
        return value > 0;
    }
    
    // GCD and LCM
    BigInt gcd(const BigInt& other) const {
        mpz_class result;
        mpz_gcd(result.get_mpz_t(), value.get_mpz_t(), other.value.get_mpz_t());
        return BigInt(result);
    }
    
    BigInt lcm(const BigInt& other) const {
        mpz_class result;
        mpz_lcm(result.get_mpz_t(), value.get_mpz_t(), other.value.get_mpz_t());
        return BigInt(result);
    }
    
    // Modular inverse
    BigInt mod_inverse(const BigInt& modulus) const {
        mpz_class result;
        if (mpz_invert(result.get_mpz_t(), value.get_mpz_t(), modulus.value.get_mpz_t()) == 0) {
            throw std::runtime_error("Modular inverse does not exist");
        }
        return BigInt(result);
    }
    
    // Square root
    BigInt sqrt() const {
        if (value < 0) {
            throw std::runtime_error("Square root of negative number");
        }
        mpz_class result;
        mpz_sqrt(result.get_mpz_t(), value.get_mpz_t());
        return BigInt(result);
    }
    
    // Access to underlying mpz_class for advanced operations
    const mpz_class& get_mpz() const {
        return value;
    }
    
    mpz_class& get_mpz() {
        return value;
    }
    
    // Random number generation
    static BigInt random(const BigInt& max) {
        static gmp_randclass rng(gmp_randinit_default);
        static bool seeded = false;
        if (!seeded) {
            rng.seed(time(nullptr));
            seeded = true;
        }
        return BigInt(rng.get_z_range(max.value));
    }
    
    // Stream operators
    friend std::ostream& operator<<(std::ostream& os, const BigInt& bi) {
        os << bi.value;
        return os;
    }
    
    friend std::istream& operator>>(std::istream& is, BigInt& bi) {
        is >> bi.value;
        return is;
    }
};

// Global operators for mixed types
inline BigInt operator+(long lhs, const BigInt& rhs) {
    return BigInt(lhs) + rhs;
}

inline BigInt operator-(long lhs, const BigInt& rhs) {
    return BigInt(lhs) - rhs;
}

inline BigInt operator*(long lhs, const BigInt& rhs) {
    return BigInt(lhs) * rhs;
}

inline BigInt operator/(long lhs, const BigInt& rhs) {
    return BigInt(lhs) / rhs;
}

inline BigInt operator%(long lhs, const BigInt& rhs) {
    return BigInt(lhs) % rhs;
}

inline bool operator<(long lhs, const BigInt& rhs) {
    return BigInt(lhs) < rhs;
}

inline bool operator>(long lhs, const BigInt& rhs) {
    return BigInt(lhs) > rhs;
}

inline bool operator<=(long lhs, const BigInt& rhs) {
    return BigInt(lhs) <= rhs;
}

inline bool operator>=(long lhs, const BigInt& rhs) {
    return BigInt(lhs) >= rhs;
}

inline bool operator==(long lhs, const BigInt& rhs) {
    return BigInt(lhs) == rhs;
}

inline bool operator!=(long lhs, const BigInt& rhs) {
    return BigInt(lhs) != rhs;
}

#endif // GMP_WRAPPER_HPP
