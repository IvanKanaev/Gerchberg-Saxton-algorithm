#pragma once
#include "types.h"

#ifndef FFT_HPP
#define FFT_HPP

void swap(std::vector<std::complex<double>>& vctr, unsigned int i, unsigned int j);

void bit_reverse(std::vector<std::complex<double>>& x);

std::vector<std::complex<double>> fftfreq(const int n, const float d = 1);

std::vector<std::complex<double>> ifftfreq(const int n, const float d = 1);

std::vector<std::complex<double>> fft(const std::vector<std::complex<double>>& x);

std::vector<std::complex<double>> ifft(const std::vector<std::complex<double>>& x);

/* class fft2D {
public:
    // Прямое 2D FFT: пространственная область -> частотная
    Matrix_2DC forward(const Matrix_2DC& input);

    // Обратное 2D FFT: частотная область -> пространственная
    Matrix_2DC inverse(const Matrix_2DC& input);
}; */

#endif