#pragma once
#include "types.h"

void swap(std::vector<std::complex<double>>& vctr, unsigned int i, unsigned int j);

void bit_reverse(std::vector<std::complex<double>>& x);

std::vector<std::complex<double>> fftfreq(const int n);

std::vector<std::complex<double>> ifftfreq(const int n);

std::vector<std::complex<double>> fft(const std::vector<std::complex<double>>& x);

std::vector<std::complex<double>> ifft(const std::vector<std::complex<double>>& x);

Matrix_2DC fft2d(const Matrix_2DC& x);

Matrix_2DC ifft2d(const Matrix_2DC& x);