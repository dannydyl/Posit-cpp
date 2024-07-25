# Posit Arithmetic Library in C++

## Overview

This C++ library is specifically designed to facilitate the implementation and evaluation of the posit number system, a new alternative to floating-point arithmetic that offers improved accuracy, performance, and reproducibility across different computing environments. It supports multiple configurations with posit sizes of 8, 16, and 32 bits and allows for any exponent size (es), catering to various precision requirements. The library is suitable for educational purposes, simulations, and research in numerical methods, digital signal processing, and high-performance computing.

## What is a Posit?

Posit is a type of number designed to be a direct replacement for IEEE standard floating points. Unlike traditional floating points, posits provide tapered precision when numbers are close to zero or very large, which can lead to more accurate and robust calculations, particularly in applications requiring high dynamic range or where rounding errors are problematic. Posits use a combination of a sign bit, regime bits, exponent bits, and fraction bits to represent real numbers, adjusting the number of bits used for precision based on the magnitude of the value.

## Files Description

### `posit_util.h`
This header file defines the framework and necessary functions for operating with posits. It encapsulates the following features:
- Templates for creating posits with customizable sizes and exponent bits.
- Functions for basic arithmetic operations: addition, subtraction, multiplication, and division.
- Comparative operations to facilitate sorting and algorithmic decisions.

### `posit_util.cpp`
Implements the computational logic described in `posit_util.h`. This source file is crucial for executing the arithmetic operations, ensuring that all posit calculations are performed according to the specified configurations and accurately reflecting the unique properties of the posit number system.

### `main.cpp`
Provides examples on initializing posits and performing arithmetic operations. It serves as a practical guide to using the library and can be used to validate the implementation by comparing results with theoretical expectations or other arithmetic systems.
