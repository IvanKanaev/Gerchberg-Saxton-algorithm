#pragma once
#include "types.h"

class FFT2D {
public:
    // Прямое 2D FFT: пространственная область -> частотная
    Matrix_2DC forward(const Matrix_2DC& input);

    // Обратное 2D FFT: частотная область -> пространственная
    Matrix_2DC inverse(const Matrix_2DC& input);
};