#pragma once
#include "types.h"
#include "fft.h"

/**
 * Вычисляет передаточную функцию H(fx, fy) = exp(-i * k * z * sqrt(1 - (lambda*fx)^2 - (lambda*fy)^2)).
 *
 * @param nx, ny        Размер сетки (строки, столбцы).
 * @param dx, dy        Шаг пикселя в пространственной области (метры).
 * @param wavelength     Длина волны (метры).
 * @param z              Расстояние распространения (метры).
 * @return               Матрица комплексных чисел H размером nx × ny.
 */
Matrix_2DC compute_transfer_function(
    size_t nx, size_t ny,
    double dx, double dy,
    double wavelength,
    double z
);

/**
 * Распространяет комплексное поле `field` на расстояние `z`.
 * Использует переданный объект FFT2D для прямого и обратного преобразований.
 *
 * @param field          Исходное комплексное поле.
 * @param dx, dy         Шаг пикселя.
 * @param wavelength     Длина волны.
 * @param z              Расстояние распространения.
 * @param fft            Экземпляр FFT2D (по ссылке, чтобы использовать один объект).
 * @return               Поле после распространения.
 */
Matrix_2DC propagate(
    const Matrix_2DC& field,
    double dx, double dy,
    double wavelength,
    double z,
    FFT2D& fft
);