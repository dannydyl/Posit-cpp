#include <iostream>
#include <bitset>
#include "posit_util.h"

int main() {
    // Example POSIT values using bitsets
    std::bitset<8> posit1Bitset("01010000"); // POSIT1:
    std::bitset<8> posit2Bitset("10110000"); // POSIT2:

    // Convert bitsets to vectors
    std::vector<bool> posit1 = MyPosit::bitsetToVector8(posit1Bitset);
    std::vector<bool> posit2 = MyPosit::bitsetToVector8(posit2Bitset);

    int es = 1; // Exponent size

    MyPosit p1(posit1, 8, es);
    MyPosit p2(posit2, 8, es);

    std::cout << "-----------------p1 double representation: " << std::endl;
    p1.printDouble();

    std::cout << "-----------------p2 double representation: " << std::endl;
    p2.printDouble();

    // Perform addition of p1 and p2
    MyPosit result = p1.add(p1, p2);
    std::cout << "Result of p1 + p2:" << std::endl;
    result.printDouble();

    // MyPosit result = p1.mul(p1, p2);
    // std::cout << "Result of p1 * p2:" << std::endl;
    // result.printDouble();

    // how to run
    // g++ -std=c++20 main.cpp posit_util.cpp -o main
    // ./main
    return 0;
}