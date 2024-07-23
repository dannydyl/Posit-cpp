#ifndef POSIT_UTIL_H
#define POSIT_UTIL_H

#include <vector>
#include <bitset>

class MyPosit{
public:
    MyPosit(const std::vector<bool>& posit, int nbits, int es); 
    MyPosit add(MyPosit& p1, MyPosit& p2);
    MyPosit sub(MyPosit& p1, MyPosit& p2);
    MyPosit mul(const MyPosit& p1, const MyPosit& p2);
    MyPosit div(const MyPosit& p1, const MyPosit& p2);
    void printComponents() const;
    void printDouble();
    double getDouble();
    static std::vector<bool> bitsetToVector8(const std::bitset<8>& bitset);
    static std::vector<bool> bitsetToVector16(const std::bitset<16>& bitset);
    static std::vector<bool> bitsetToVector32(const std::bitset<32>& bitset);

private:
    std::vector<bool> positValue;
    bool sign;
    int regimeValue;
    std::vector<bool> regimeBits;   // no need
    std::vector<bool> exponent;
    std::vector<bool> fraction; // no need
    uint32_t fractionInt;
    int es; // exponent size
    int scalingFactor;
    int useed;  // no need
    int32_t effectiveExponent;
    int absoluteEffectiveExponent;  // no need
    int nbits;
    bool isZeroFlag;
    bool isNaRFlag;

    static int countLeadingBits(const std::vector<bool>& bits, int nbits);
    void extractPosit(const std::vector<bool>& posit, int nbits, int es);
    void calculateScalingFactor();
    void calculateEffectiveExponent();
    void calculateAbsoluteEffectiveExponent();
    void normalizeResultFraction(MyPosit& result);
    int clog2(int width);

    // Bit shift functions
    uint32_t leftShift(uint32_t value, int shift);
    uint32_t rightShift(uint32_t value, int shift);
    uint64_t leftShift(uint64_t value, int shift);
    uint64_t rightShift(uint64_t value, int shift);

    // checking NaR, Zero
    void isNaR(const std::vector<bool>& posit, int nbits);
    void isZero(const std::vector<bool>& posit, int nbits);
    MyPosit getNaR(int nbits, int es);
    MyPosit getZero(int nbits, int es);

    // hidden bit
    uint32_t hiddenBit(uint32_t fraction);

    // count leading zeros
    int countLeadingZeros(int n);

    // floor division
    int floorDivision(int a, int b);

    // masking function
    uint32_t mask(int bits, int size);

    // two's complement function
    uint32_t twosComplement(uint32_t value, int nbits);

    // operations (scaling embedded)
    MyPosit add_internal(const MyPosit& p1, const MyPosit& p2);
    MyPosit sub_internal(const MyPosit& p1, const MyPosit& p2);
    MyPosit mul_internal(const MyPosit& p1, const MyPosit& p2);
    MyPosit div_internal(const MyPosit& p1, const MyPosit& p2);

    // decode
    std::vector<bool> encodePosit(MyPosit decodedPosit, int nbits, int es);
    std::vector<bool> toVectorBool(uint32_t encodedPosit, int nbits);
};

#endif // POSIT_UTIL_H
