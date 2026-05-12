#pragma once
#include "types.h"

/**
 * Дополняет матрицу нулями до квадратной.
 * Размер стороны квадрата = max(rows, cols).
 */
Matrix_2D pad_to_square(const Matrix_2D& mat);

/**
 * Изменяет размер матрицы до new_size × new_size.
 * Использует билинейную интерполяцию.
 */
Matrix_2D resize_matrix(const Matrix_2D& mat, unsigned int new_size);

