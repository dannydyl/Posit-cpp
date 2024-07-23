#include <iostream>
#include <bitset>
#include <cmath>
#include <vector>
#include <climits>
#include <algorithm>
#include "posit_util.h"

MyPosit::MyPosit(const std::vector<bool>& posit, int nbits, int es) : positValue(posit), nbits(nbits), es(es){
    extractPosit(posit, nbits, es);
    calculateScalingFactor();
    calculateEffectiveExponent();
}

// Print components for debugging
void MyPosit::printComponents() const {
    std::cout << "Sign: " << sign << std::endl;
    std::cout << "regime bits: " << std::endl; 
    for (bool bit : regimeBits) std::cout << bit;
    std::cout << std::endl;
    std::cout << "regimeValue : " << regimeValue << std::endl;
    std::cout << "Exponent: ";
    for (bool bit : exponent) std::cout << bit;
    std::cout << std::endl;
    std::cout << "Fraction: ";
    for (bool bit : fraction) std::cout << bit;
    std::cout << std::endl;
    std::cout << "Fraction Int: " << fractionInt << std::endl;
    std::cout << "Scaling Factor: " << scalingFactor << std::endl; // Debug output
    std::cout << "Effective Exponent: " << effectiveExponent << std::endl;
}

double MyPosit::getDouble(){

    // std::cout << "\n Get double function \n" << std::endl;
    if (isNaRFlag) {
        // std::cout << "NaR detected" << std::endl;
        return 0.0;
    } else if (isZeroFlag) {
        // std::cout << "Zero detected" << std::endl;
        return 0.0;
    }

    int fexp = effectiveExponent + 1023; // floating exponent

    int width = 32;
    uint64_t fexpbits;
    uint64_t ffracbits;
    uint32_t msb = 1U << 7;

    if (fexp > 2046){ // overflow
        fexpbits = leftShift(2046U, 53);
        ffracbits = std::numeric_limits<uint64_t>::max();
        // std::cout << "overflow" << std::endl;
    } else if (fexp < 1){   // underflow
        fexpbits = 0;
        if (width <= 64) {
            ffracbits = leftShift(static_cast<uint64_t>(msb | rightShift(fractionInt, 1)), 64 - width);
        } else {
            ffracbits = rightShift(msb | rightShift(fractionInt, 1), width - 64);
        }
        ffracbits = rightShift(ffracbits, -fexp);
        // std::cout << "underflow" << std::endl;
    } else {    // normal
        fexpbits = leftShift(static_cast<uint64_t>(fexp & 0x7FF), 53);
        if (width <= 32){
            ffracbits = leftShift(static_cast<uint64_t>(fractionInt), 64 - width);
        } else {
            ffracbits = rightShift(fractionInt, width - 64);
        }
        // std::cout << "normal" << std::endl;
        // std::cout << "ffracbits (normal): " << ffracbits << std::endl;
    }
    
    uint64_t result = ffracbits;
    result = fexpbits | rightShift(result, 11);
    result = leftShift(static_cast<uint64_t>(sign), 63) | rightShift(result, 1);
    
    // don't underflow to zero
    if (leftShift(result, 1) == 0) {
        result++;
    }
    
    // Convert the uint64_t to a double
    double finalValue ;
    std::memcpy(&finalValue, &result, sizeof(double));
    
    // Print the result
    // std::cout << "fraction Int : " << fractionInt << std::endl;
    // std::cout << "fexp : " << fexp << std::endl;
    // std::cout << "ffracbits : " << ffracbits << std::endl;
    // std::cout << "result : " << result << std::endl;
    // std::cout << "fexpbits : " << fexpbits << std::endl;
    // std::cout << "finalValue : " << finalValue << std::endl;

    return finalValue;
}

std::vector<bool> MyPosit::encodePosit(MyPosit decodedPosit, int nbits, int es){
    uint32_t encode;
    uint32_t regimebits;
    uint32_t exponentbits;

    int maxexp = (1 << es) * (nbits - 2);   // maximum exponent is 2^(es) * (nbits - 2)
    if(decodedPosit.effectiveExponent < -maxexp){   // hanndling underflow
        decodedPosit.effectiveExponent = -maxexp;
    } else if(decodedPosit.effectiveExponent > maxexp){ // handling overflow
        decodedPosit.effectiveExponent = maxexp;
    }
    // std::cout << "decodedPosit.effectiveExponent: " << decodedPosit.effectiveExponent << std::endl;
    int regimeInt = floorDivision(decodedPosit.effectiveExponent, (1 << es)); // regime value
    // std::cout << "regimebits declaration: " << regimebits << std::endl;
    int regimeLength = std::max(-regimeInt + 1, regimeInt + 2); // regime length

    if (1 + regimeLength + es >= nbits && decodedPosit.fractionInt >= 0x80000000){
        decodedPosit.effectiveExponent++;
        regimeInt = floorDivision(decodedPosit.effectiveExponent, (1 << es));
        regimeLength = std::max(-regimeInt + 1, regimeInt + 2);
    }

    uint32_t encodedEffectiveExponent = decodedPosit.effectiveExponent - (1 << es) * regimeInt;

    if(regimeInt < 0){
        regimebits = rightShift(0x80000000, -regimeInt);
    } else {
        regimebits =  mask(0xFFFFFFFF, regimeInt + 1); 
    }

    exponentbits = mask(leftShift(encodedEffectiveExponent, 32 - es), es);

    encode = decodedPosit.fractionInt;
    encode = exponentbits | rightShift(encode, es);
    encode = regimebits | rightShift(encode, regimeLength);
    encode = rightShift(encode, 1);
    // std::cout << "before encode int value: " << encode << std::endl;
    if(decodedPosit.sign){
        encode = twosComplement(encode, nbits);
    } else {
        encode = mask(encode, nbits);
    }
    // std::cout << "after encode int value: " << encode << std::endl;
    // std::cout << "regimeInt: " << regimeInt << std::endl;
    // std::cout << "regimeLength: " << regimeLength << std::endl;
    // std::cout << "encode regimebits: " << regimebits << std::endl;
    // std::cout << "encode exponentbits: " << exponentbits << std::endl;
    // std::cout << "encode int value: " << encode << std::endl;
    return toVectorBool(encode, nbits);
}   

std::vector<bool> MyPosit::toVectorBool(uint32_t value, int nbits){
    std::vector<bool> result(nbits, false);
    int startShift = 32 - nbits; // Calculate the shift to start from MSB
    for (int i = 0; i < nbits; ++i) {
        // Shift the value right by 'nbits-i-1' positions and extract the least significant bit
        result[i] = (value >> (startShift + nbits - i - 1)) & 1;
    }
    // std::cout << "@@@@@@@@@@@Result vector: ";
    // for (bool bit : result) {
    //     std::cout << bit;
    // }
    // std::cout << std::endl;
    return result;
}

void MyPosit::printDouble() {
    double temp = getDouble();
    std::cout << "Value in double type : " << temp << std::endl;
}

// Function to count leading zeros or ones based on the first bit after the sign bit
int MyPosit::countLeadingBits(const std::vector<bool>& bits, int nbits) {
    bool countLeadingBits = bits[0];  // Set the value to count based on bit at index 6
    int count = 0;

    for (int i = 0; i <= nbits; ++i) {
        if (bits[i] == countLeadingBits) {
            ++count;
        } else {
            break;
        }
    }
    
    return count;
}

// Function to extract POSIT components
void MyPosit::extractPosit(const std::vector<bool>& posit, int nbits, int es) {
    this->sign = posit[0];
    std::vector<bool> unsignedPosit(nbits - 1);  // Store the unsigned POSIT value
    for (size_t i = 1; i < nbits; ++i) { 
        unsignedPosit[i-1] = posit[i];
    }    

    // Check for NaR and Zero
    isNaR(posit, nbits);
    isZero(posit, nbits);

    // std::cout << "111 from extractPosit unsignedPosit: ";
    // for (bool bit : unsignedPosit) {
    //     std::cout << bit;
    // }
    // std::cout << std::endl;

    if (this->sign == 1) {
        // Perform two's complement for 7-bit unsignedPosit if sign bit is 1
        // Flip all bits for two's complement
        // std::cout << "negative posit" << std::endl;
        for (size_t i = 0; i < unsignedPosit.size(); ++i) {
             unsignedPosit[i] = !unsignedPosit[i];
        }
        // unsignedPosit = (unsignedPosit + 1) & 0x7F;  // Add 1 and ensure it stays within 7 bits
        bool carry = true; // Start by adding 1
        for (size_t i = unsignedPosit.size() - 1; i >= 0; --i) {
            bool newBit = unsignedPosit[i] ^ carry; // XOR current bit with carry
            carry = unsignedPosit[i] && carry;      // AND current bit with carry (propagate the carry if both are 1)
            unsignedPosit[i] = newBit;
            if (!carry) break; // If no carry, no need to continue
        }

        if (unsignedPosit.size() > nbits - 1) {
            unsignedPosit.resize(nbits - 1); 
        }

    }

    // std::cout << "from extractPosit unsignedPosit after twos complement: ";
    // for (bool bit : unsignedPosit) {
    //     std::cout << bit;
    // }
    // std::cout << std::endl;

    // Determine if we are counting leading zeros or ones
    int leadingBits = countLeadingBits(unsignedPosit, nbits);
    // std::cout << "-------Leading Bits in extract posit: " << leadingBits << std::endl;

    // Determine regimeValue
    this->regimeValue = (unsignedPosit[0] == 0) ? -leadingBits : leadingBits - 1;
    // std::cout << "-------Regime Value in extract posit: " << this->regimeValue << std::endl;

    // Extract exponent
    int exponentStart = leadingBits + 1; // Move to the bit after the terminating bit
    // std::vector<bool> exponent(es, 0);
    this->exponent.resize(es, 0);
    for (size_t i = 0; i < es ; ++i) {
        this->exponent[i] = unsignedPosit[exponentStart - i];
    }
    // std::cout << "Exponent start: " << exponentStart << std::endl;
    // for(bool bit : exponent){
    //     std::cout << bit;
    // }
    // Extract fraction
    int fractionStart = exponentStart + es;
    this->fraction.clear();
    // std::cout << "Fraction start: " << fractionStart << std::endl;

    for (int i = fractionStart; i <= unsignedPosit.size(); ++i) {
        fraction.push_back(unsignedPosit[i]);
    }

    // Calculate fractionInt from fraction bits
    this->fractionInt = 0;
    for (int i = fraction.size() - 1; i >= 0; --i) {
        if (fraction[i]) {
            this->fractionInt |= (1U << (fraction.size() - 1 - i));
        }
    }
    this->fractionInt = leftShift(fractionInt, 32 - fraction.size());   // left shift to the point where the fraction starts in total 32 bit integer
    // std::cout << "fraction int from extractPosit: " << fractionInt << std::endl;
}

void MyPosit::calculateScalingFactor() {
    this->useed = std::pow(2 , std::pow(2, es));
    this->scalingFactor = std::pow(useed, regimeValue);
    // std::cout << "Scaling factor: " << scalingFactor << std::endl;
}

void MyPosit::calculateEffectiveExponent(){
    this->effectiveExponent = this->regimeValue * (1 << es);
    for (size_t i = 0; i < exponent.size(); ++i) {
        this->effectiveExponent += exponent[i] * (1 << (es - i - 1));
    }
    // std::cout << "Effective Exponent: " << effectiveExponent << std::endl;
    // std::cout << "regimeValue: " << regimeValue << std::endl;
}

int MyPosit::clog2(int width) {
    return static_cast<int>(std::ceil(std::log2(width)));
}

void MyPosit::calculateAbsoluteEffectiveExponent(){
    int r_width = clog2(8); // assuming POSIT8 format
    int msb = (this->effectiveExponent >> (this->es + r_width)) & 1;
    int mask = msb ? ~((1 << (es + r_width + 1)) - 1) : 0;
    this->absoluteEffectiveExponent = (mask ^ this->effectiveExponent) + msb;
}


MyPosit MyPosit::add(MyPosit& p1, MyPosit& p2) {
    if (p1.isNaRFlag || p2.isNaRFlag) {
        return getNaR(p1.nbits, p1.es);
    } else if (p1.isZeroFlag || p2.isZeroFlag) {
        // Return p2 if p1 is zero
        if (p1.isZeroFlag) {
            return p2;  // Ensure p2 is returned exactly as it is
        }
        // Return p1 if p2 is zero
        else {
            return p1;  // Ensure p1 is returned exactly as it is
        }
    } else if (p1.sign == p2.sign){
        return add_internal(p1, p2);
    } else if (p1.getDouble() == -p2.getDouble()){
        return getZero(p1.nbits, p1.es);
    } else {
        // subtract
        return sub_internal(p1, p2);
    }
}

MyPosit MyPosit::sub(MyPosit& p1, MyPosit& p2) {
    if (p1.isNaRFlag || p2.isNaRFlag) {
        return getNaR(p1.nbits, p1.es);
    } else if (p1.isZeroFlag || p2.isZeroFlag) {
        // Return p2 if p1 is zero
        if (p1.isZeroFlag) {
            return p2;  // Ensure p2 is returned exactly as it is
        }
        // Return p1 if p2 is zero
        else {
            return p1;  // Ensure p1 is returned exactly as it is
        }
    } else if(p1.sign == p2.sign){
        return sub_internal(p1, p2);
    } else if (p1.getDouble() == p2.getDouble()){
        return getZero(p1.nbits, p1.es);
    } else {
        // subtract
        return add_internal(p1, p2);
    }
}

MyPosit MyPosit::mul(const MyPosit& p1, const MyPosit& p2) {
    if(p1.isNaRFlag || p2.isNaRFlag){
        return getNaR(p1.nbits, p1.es);
    } else if (p1.isZeroFlag || p2.isZeroFlag){
        return getZero(p1.nbits, p1.es);
    } else {
        return mul_internal(p1, p2);
    }
}

MyPosit MyPosit::div(const MyPosit& p1, const MyPosit& p2) {
    if(p1.isNaRFlag || p2.isNaRFlag){
        return getNaR(p1.nbits, p1.es);
    } else if (p1.isZeroFlag || p2.isZeroFlag){
        return getZero(p1.nbits, p1.es);
    } else {
        return div_internal(p1, p2);
    }
}

// Add two POSIT numbers
MyPosit MyPosit::add_internal(const MyPosit& p1, const MyPosit& p2) {
    // hidden bit
    uint64_t p1FractionInt = hiddenBit(p1.fractionInt);
    uint64_t p2FractionInt = hiddenBit(p2.fractionInt);
    uint64_t resultFractionInt;
    // std::cout << "p1Copy.fractionInt after hidden bit: " << p1FractionInt << std::endl;
    // std::cout << "p2Copy.fractionInt after hidden bit: " << p2FractionInt << std::endl;
    MyPosit result(p1.positValue, p1.nbits, p1.es);

    // Shift the fraction of the posit with the smaller exponent
    if (p1.effectiveExponent > p2.effectiveExponent) {
        // std::cout << "p1 is greater" << std::endl;
        result.effectiveExponent = p1.effectiveExponent;
        // Shift p2Copy's fraction
        p2FractionInt = rightShift(p2FractionInt, p1.effectiveExponent - p2.effectiveExponent);
    } else {
        // std::cout << "p2 is greater" << std::endl;
        result.effectiveExponent = p2.effectiveExponent;
        // Shift p1Copy's fraction
        p1FractionInt = rightShift(p1FractionInt, p2.effectiveExponent - p1.effectiveExponent);
    }


    resultFractionInt = p1FractionInt + p2FractionInt;
    // std::cout << "p1Copy.fractionInt after shift: " << p1FractionInt << std::endl;
    // std::cout << "p2Copy.fractionInt after shift: " << p2FractionInt << std::endl;
    // std::cout << "result fraction int: " << resultFractionInt << std::endl;
    if (rightShift(resultFractionInt, 32) != 0) {
        // std::cout << "Overflow detected" << std::endl;
        result.effectiveExponent++;
        resultFractionInt = rightShift(resultFractionInt, 1);
    }
    
    result.sign = p1.sign;
    result.fractionInt = leftShift(resultFractionInt, 1);
    // std::cout << "result.fractionInt after overflow check: " << result.fractionInt << std::endl;
    
    result.positValue = encodePosit(result, p1.nbits, p1.es);
    MyPosit res(result.positValue, p1.nbits, p1.es);
    return res;
}

MyPosit MyPosit::sub_internal(const MyPosit& p1, const MyPosit& p2) {
    // std::cout << "@@@@@@@@@@@Subtracting" << std::endl;
    // hidden bit
    uint64_t p1FractionInt = hiddenBit(p1.fractionInt);
    uint64_t p2FractionInt = hiddenBit(p2.fractionInt);
    uint64_t resultFractionInt;
    // std::cout << "p1Copy.fractionInt after hidden bit: " << p1FractionInt << std::endl;
    // std::cout << "p2Copy.fractionInt after hidden bit: " << p2FractionInt << std::endl;

    // initialize the result using first Posit's bitset and properities
    MyPosit result(p1.positValue, p1.nbits, p1.es);

    // Determine which Posit has the greater magnitude by comparing exponents and fractions
    if(p1.effectiveExponent > p2.effectiveExponent || (p1.effectiveExponent == p2.effectiveExponent && p1.fractionInt > p2.fractionInt)){
        // std::cout << "p1 is greater" << std::endl;
        result.sign = p1.sign;
        result.effectiveExponent = p1.effectiveExponent;
        // Shift p2Copy's fraction
        p2FractionInt = rightShift(p2FractionInt, p1.effectiveExponent - p2.effectiveExponent);
        resultFractionInt = p1FractionInt - p2FractionInt;
    } else {
        // std::cout << "p2 is greater" << std::endl;
        result.sign = !p1.sign;
        result.effectiveExponent = p2.effectiveExponent;
        // Shift p1Copy's fraction
        p1FractionInt = rightShift(p1FractionInt, p2.effectiveExponent - p1.effectiveExponent);
        resultFractionInt = p2FractionInt - p1FractionInt;
    }

    result.effectiveExponent -= countLeadingZeros(resultFractionInt);
    result.fractionInt = leftShift(resultFractionInt, countLeadingZeros(resultFractionInt) + 1);
    result.positValue = encodePosit(result, p1.nbits, p1.es);
    // for(bool bit : result.positValue){
    //     std::cout << bit;
    // }
    // std::cout << "@@@" << std::endl;
    MyPosit res(result.positValue, p1.nbits, p1.es);
    return res;
}

MyPosit MyPosit::mul_internal(const MyPosit& p1, const MyPosit& p2) {
    uint64_t p1FractionInt = hiddenBit(p1.fractionInt);
    uint64_t p2FractionInt = hiddenBit(p2.fractionInt);
    // std::cout << "@@@@@@@@@@@@@multiplying" << std::endl;
    // initialize the result using first Posit's bitset and properities
    MyPosit result(p1.positValue, p1.nbits, p1.es);

    // Multiply the two fractions with hidden bits and normalize by shifting right by 32.
    // This adjusts the result's magnitude back to the appropriate range.
    uint64_t resultFractionInt = rightShift(p1FractionInt * p2FractionInt, 32);
    result.effectiveExponent = p1.effectiveExponent + p2.effectiveExponent + 1;

    // If the hidden bit is not present (resultFractionInt is not normalized), decrement the exponent
    // and shift the fraction to the left to restore the normalization.
    if ((resultFractionInt & 0x80000000) == 0){
        result.effectiveExponent--;
        resultFractionInt = leftShift(resultFractionInt, 1);
    }

    // Determine the sign of the result based on the XOR of the signs of the inputs.
    result.sign = p1.sign ^ p2.sign;
    result.fractionInt = leftShift(resultFractionInt, 1);

    result.positValue = encodePosit(result, p1.nbits, p1.es);
    MyPosit res(result.positValue, p1.nbits, p1.es);
    return res;
}

MyPosit MyPosit::div_internal(const MyPosit& p1, const MyPosit& p2) {
    uint64_t p1FractionInt = hiddenBit(p1.fractionInt);
    uint64_t p2FractionInt = hiddenBit(p2.fractionInt);
    MyPosit result(p1.positValue, p1.nbits, p1.es);

    result.effectiveExponent = p1.effectiveExponent - p2.effectiveExponent;

    if (p1FractionInt < p2FractionInt){
        result.effectiveExponent--;
        p2FractionInt = rightShift(p2FractionInt, 1);
    }

    result.sign = p1.sign ^ p2.sign;
    result.fractionInt = leftShift(p1FractionInt, 32) / p2FractionInt;

    result.positValue = encodePosit(result, p1.nbits, p1.es);
    MyPosit res(result.positValue, p1.nbits, p1.es);
    return res;
}

// Overloaded leftShift function
uint32_t MyPosit::leftShift(uint32_t value, int shift) {
    if (shift >= 0) {
        return value << shift;
    } else {
        return 0;
    }
}

uint64_t MyPosit::leftShift(uint64_t value, int shift) {
    if (shift >= 0) {
        return value << shift;
    } else {
        return 0;
    }
}

// Overloaded rightShift function
uint32_t MyPosit::rightShift(uint32_t value, int shift) {
    if (shift >= 0) {
        return value >> shift;
    } else {
        return 0;
    }
}

uint64_t MyPosit::rightShift(uint64_t value, int shift) {
    if (shift >= 0) {
        return value >> shift;
    } else {
        return 0;
    }
}

uint32_t MyPosit::hiddenBit(uint32_t fractionInt) {
    return rightShift(fractionInt, 1) | 0x80000000;
}

int MyPosit::countLeadingZeros(int n) {
    if (n == 0) {
        return CHAR_BIT * sizeof(n);  // CHAR_BIT is the number of bits per byte, usually 8
    } else {
        return __builtin_clz(n);  // Use GCC's built-in function to count leading zeros
    }
}

int MyPosit::floorDivision(int a, int b) {
    return ((a / b) - ((a) % (b) < 0));
}

void MyPosit::isNaR(const std::vector<bool>& posit, int nbits) {
    // Create a NaR representation for the given number of bits
    std::vector<bool> nar(nbits, false); // Start with all bits as false
    nar[0] = true; // Set the first bit to true (NaR)

    // Directly compare the given posit vector with the NaR vector
    this->isNaRFlag = (posit == nar);
    if(this->isNaRFlag){
        // std::cout << "NaR detected" << std::endl;
    }
}

void MyPosit::isZero(const std::vector<bool>& posit, int nbits) {
    // Create a zero representation for the given number of bits
    std::vector<bool> zero(nbits, false); // Start with all bits as false
    
    // Directly compare the given posit vector with the zero vector
    this->isZeroFlag = (posit == zero);
    if(this->isZeroFlag){
        // std::cout << "Zero detected" << std::endl;
    }
}

MyPosit MyPosit::getNaR(int nbits, int es) {
    // Create a NaR representation for the given number of bits
    std::vector<bool> nar(nbits, false); // Start with all bits as false
    nar[0] = true; // Set the first bit to true (NaR)

    return MyPosit(nar, nbits, es);
}

MyPosit MyPosit::getZero(int nbits, int es) {
    // Create a zero representation for the given number of bits
    std::vector<bool> zero(nbits, false); // Start with all bits as false

    return MyPosit(zero, nbits, es);
}

uint32_t MyPosit::mask(int bits, int size) {
    return ((bits) & leftShift(0xFFFFFFFF, 32 - (size)));
}

uint32_t MyPosit::twosComplement(uint32_t value, int nbits) {
    return mask(-mask(value, nbits), nbits);
}

std::vector<bool> MyPosit::bitsetToVector8(const std::bitset<8>& bitset){
    std::vector<bool> vec(bitset.size());
    for(size_t i = 0; i < bitset.size(); i++){
        vec[bitset.size() - 1 - i] = bitset[i];
    }
    return vec;
}

std::vector<bool> MyPosit::bitsetToVector16(const std::bitset<16>& bitset){
    std::vector<bool> vec(bitset.size());
    for(size_t i = 0; i < bitset.size(); i++){
        vec[bitset.size() - 1 - i] = bitset[i];
    }
    return vec;
}

std::vector<bool> MyPosit::bitsetToVector32(const std::bitset<32>& bitset){
    std::vector<bool> vec(bitset.size());
    for(size_t i = 0; i < bitset.size(); i++){
        vec[bitset.size() - 1 - i] = bitset[i];
    }
    return vec;
}