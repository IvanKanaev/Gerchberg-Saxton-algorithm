#pragma once
#include <vector>
#include <complex>

// матрица для действительных чисел (интенсивность, амплитуда, фаза и т.д.)
using Matrix_2D = std::vector<std::vector<double>>;

// матрица для комплексных чисел 
using Matrix_2DC = std::vector<std::vector<std::complex<double>>>;