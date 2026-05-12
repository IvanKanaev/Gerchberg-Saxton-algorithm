#pragma once
#include <vector>
#include "types.h"
#include "fft.h"

/**
 * Результат работы алгоритма: амплитуда и фаза в начальной плоскости.
 */
struct GSResult {
    Matrix_2D amplitude;   // восстановленная амплитуда (нормированная)
    Matrix_2D phase;       // восстановленная фаза в радианах
};

/**
 * Алгоритм Герчберга–Сакстона (Gerchberg–Saxton).
 *
 * Восстанавливает фазу по набору измеренных интенсивностей
 * на разных расстояниях вдоль оси z.
 *
 * @param intensities  Набор квадратных матриц интенсивности (одинакового размера).
 *                     Каждая матрица соответствует одной плоскости.
 * @param z_distances  Координаты z этих плоскостей (в метрах, по возрастанию).
 * @param wavelength   Длина волны (в метрах).
 * @param dx, dy       Шаг пикселя (в метрах).
 * @param fft          Экземпляр FFT2D (реализация кента).
 * @param max_iter     Максимальное количество итераций (по умолчанию 100).
 * @param tolerance    Порог сходимости (по умолчанию 1e-6).
 * @return             Структура GSResult с amplitude и phase.
 */
GSResult gerchberg_saxton(
    const std::vector<Matrix_2D>& intensities,
    const std::vector<double>& z_distances,
    double wavelength,
    double dx, double dy,
    FFT2D& fft,
    int max_iter = 100,
    double tolerance = 1e-6
);