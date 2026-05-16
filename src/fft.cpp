#include "fft.h"
#include "types.h"
constexpr double PI = 3.14159265358979323846;

/**
 * @brief Обмен двух элементов в векторе комплексных чисел.
 * @param vctr   ссылка на вектор
 * @param i      индекс первого элемента
 * @param j      индекс второго элемента
 */
void swap(std::vector<std::complex<double>>& vctr, unsigned int i, unsigned int j){
    std::complex<double> buffer;
    buffer = vctr[i];
    vctr[i] = vctr[j];
    vctr[j] = buffer;
}

/**
 * @brief Бит-реверсивная перестановка элементов вектора.
 *        Используется для подготовки данных к БПФ с прореживанием по времени.
 * @note  Корректна только для длины, являющейся степенью двойки.
 * @param x вектор, который будет переставлен на месте (in-place)
 */
void bit_reverse(std::vector<std::complex<double>>& x){
    size_t n = x.size();
    unsigned int j = 0;
    for (unsigned int i = 0; i < n; i++){
        if (i < j) swap(x, i, j);
        unsigned int mask = n >> 1;
        while (j & mask){
            j &= ~mask;
            mask >>= 1;
        }
        j |=  mask;
    }
}

/**
 * @brief Генерация поворотных коэффициентов для прямого БПФ:
 *        w[k] = exp(-2πi * k / n)
 * @param n количество точек
 * @return вектор комплексных экспонент длины n
 */
std::vector<std::complex<double>> fftfreq(const int n){
    std::complex<double> I(0.0, 1.0); 
    std::vector<std::complex<double>> w;
    for (int i = 0; i < n; i++){
        w.push_back(std::exp(- std::complex<double>(2 * PI * i/n) * I));
    }
    return w;
}

/**
 * @brief Генерация поворотных коэффициентов для обратного БПФ:
 *        w[k] = exp(+2πi * k / n)
 * @param n количество точек
 * @return вектор комплексных экспонент длины n
 */
std::vector<std::complex<double>> ifftfreq(const int n){
    std::complex<double> I(0.0, 1.0); 
    std::vector<std::complex<double>> w;
    for (int i = 0; i < n; i++){
        w.push_back(std::exp(std::complex<double>(2 * PI * i/n) * I));
    }
    return w;
}

/**
 * @brief Одномерное прямое быстрое преобразование Фурье (БПФ).
 *        Алгоритм Кули-Тьюки с прореживанием по времени.
 *        Если длина входного вектора не степень двойки, он дополняется нулями.
 * @param x входной вектор комплексных чисел
 * @return вектор комплексного спектра такой же длины, как дополненный до степени двойки.
 *         Нормировка не производится (амплитуды растут с длиной сигнала).
 */
std::vector<std::complex<double>> fft(const std::vector<std::complex<double>>& x){
    std::vector<std::complex<double>> P = x;
    size_t n = P.size();
    unsigned int k = 1;
    int bits = 0;
    while (k < n){
        k *= 2;
        bits++;
    }
    P.resize(k, std::complex<double>(0));

    bit_reverse(P);

    int m;
    std::complex<double> P_0;
    std::complex<double> P_1;

    for (int s = 1; s < bits+1; s++){
        m = 1 << s;     // возвращает 2 в степени s
        std::vector<std::complex<double>> w = fftfreq(m);
        for (unsigned int i = 0; i < k; i+=m){
            for (int j = 0; j < m/2; j++){
                P_0 = P[i + j];
                P_1 = P[i + j + m / 2];
                P[i + j] = P_0 + w[j] * P_1;
                P[i + j + m / 2] = P_0 - w[j] * P_1;
            }
        }
    }
    return P;
}

/**
 * @brief Одномерное обратное быстрое преобразование Фурье (ОБПФ).
 *        Алгоритм аналогичен прямому, но используются сопряжённые поворотные коэффициенты,
 *        а результат делится на длину преобразования (k).
 * @param x входной вектор комплексных чисел (спектр)
 * @return восстановленный сигнал длины, равной дополненной до степени двойки.
 *         Нормировка: каждый элемент делится на k.
 */
std::vector<std::complex<double>> ifft(const std::vector<std::complex<double>>& x){
    std::vector<std::complex<double>> P = x;
    size_t n = P.size();
    unsigned int k = 1;
    int bits = 0;
    while (k < n){
        k *= 2;
        bits++;
    }
    P.resize(k, std::complex<double>(0));

    bit_reverse(P);

    int m;
    std::complex<double> P_0;
    std::complex<double> P_1;

    for (int s = 1; s < bits+1; s++){
        m = 1 << s;     // возвращает 2 в степени s
        std::vector<std::complex<double>> w = ifftfreq(m);
        for (unsigned int i = 0; i < k; i+=m){
            for (int j = 0; j < m/2; j++){
                P_0 = P[i + j];
                P_1 = P[i + j + m / 2];
                P[i + j] = P_0 + w[j] * P_1;
                P[i + j + m / 2] = P_0 - w[j] * P_1;
            }
        }
    }

    for (size_t i = 0; i < P.size(); i++)
        P[i] /= (double)k;

    return P;
}

/**
 * @brief Двумерное прямое БПФ.
 *        Вычисляется последовательным применением одномерного БПФ ко всем строкам,
 *        а затем ко всем столбцам полученной матрицы. Строки и столбцы, длина которых
 *        не является степенью двойки, автоматически дополняются нулями.
 * @param x двумерная матрица (вектор векторов) комплексных чисел
 * @return матрица спектра размером H_new × W_new, где H_new, W_new – ближайшие
 *         степени двойки для исходных высоты и ширины.
 */
Matrix_2DC fft2d(const Matrix_2DC& x){
    Matrix_2DC row_fft(x.size());
    for (size_t i = 0; i < x.size(); i++){
        row_fft[i] = fft(x[i]);
    }
    Matrix_2DC col_fft(row_fft[0].size());
    for (size_t j = 0; j < row_fft[0].size(); j++){
        for (size_t i = 0; i < row_fft.size(); i++){
            col_fft[j].push_back(row_fft[i][j]);
        }
        col_fft[j] = fft(col_fft[j]);
    }
    Matrix_2DC result(col_fft[0].size());
    for (size_t i = 0; i < col_fft[0].size(); i++){
        for (size_t j = 0; j < col_fft.size(); j++){
            result[i].push_back(col_fft[j][i]);
        }
    }

    return result;
}

/**
 * @brief Двумерное обратное БПФ.
 *        Построено аналогично прямому, но используются одномерные ОБПФ.
 *        Нормировка: каждая строка делится на W_new, каждый столбец на H_new,
 *        итоговый масштаб = 1/(W_new * H_new).
 * @param x двумерный спектр
 * @return восстановленная матрица размером H_new × W_new (дополненная нулями).
 */
Matrix_2DC ifft2d(const Matrix_2DC& x){
    Matrix_2DC row_fft(x.size());
    for (size_t i = 0; i < x.size(); i++){
        row_fft[i] = ifft(x[i]);
    }
    Matrix_2DC col_fft(row_fft[0].size());
    for (size_t j = 0; j < row_fft[0].size(); j++){
        for (size_t i = 0; i < row_fft.size(); i++){
            col_fft[j].push_back(row_fft[i][j]);
        }
        col_fft[j] = ifft(col_fft[j]);
    }
    Matrix_2DC result(col_fft[0].size());
    for (size_t i = 0; i < col_fft[0].size(); i++){
        for (size_t j = 0; j < col_fft.size(); j++){
            result[i].push_back(col_fft[j][i]);
        }
    }

    return result;
}