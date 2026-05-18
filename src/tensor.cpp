#include "tensor.h"
#include "fft.h"
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <complex>

static std::vector<std::vector<double>> read_cvs(const std::string& path){
    std::ifstream file(path);
    std::vector<std::vector<double>> matrix;
    std::string line;
    while (std::getline(file, line)) {
        std::vector<double> row;
        std::stringstream ss(line);
        std::string cell;

        while (std::getline(ss, cell, ',')) {
            row.push_back(std::stod(cell));
        }
        matrix.push_back(row);
    }
    return matrix;
}
static std::vector<std::vector<double>> resize_do2 (const std::vector<std::vector<double>>& matr){
    size_t row = matr.size();
    size_t col = matr[0].size();
    size_t max_size = std::max(row, col);

    std::vector<std::vector<double>> square(max_size, std::vector<double>(max_size, 0.0));

    size_t row_off = (max_size - row) / 2;
    size_t col_off = (max_size - col) / 2;
    for (size_t i = 0; i < row; ++i)
        for (size_t j = 0; j < col; ++j)
            square[row_off + i][col_off + j] = matr[i][j];

    size_t target = 1;
    while(target < max_size){
        target <<= 1;
    };
    if (target == max_size)
        return square;
    
    std::vector<std::vector<double>> result(target, std::vector<double>(target, 0.0));
    size_t off = (target - max_size) / 2;
    for (size_t i = 0; i < max_size; ++i)
        for (size_t j = 0; j < max_size; ++j)
            result[off + i][off + j] = square[i][j];
    return result;
}

static std::vector<double> fftfreq(size_t n, double d){
    std::vector<double> f(n);
    double df = 1.0 / (n * d);
    int half = static_cast<int>(n / 2);
    for (int i = 0; i <= half; ++i)
        f[i] = i * df;
    for (int i = half + 1; i < static_cast<int>(n); ++i)
        f[i] = (i - static_cast<int>(n)) * df;
    return f;
}

Tensor::Tensor(const std::vector<std::string>& paths, const std::vector<double> z_distance, double dx, double dy, double lam)
: z_distance(z_distance), lam(lam), dx(dx), dy(dy)
{
    for(std::vector<std::string>::const_iterator it = paths.begin(); it != paths.end(); it++){
        data.push_back(resize_do2(read_cvs(*it)));
    };
    this->nx = data[0].size();
    this->ny = data[0][0].size();
    this->fx = fftfreq(nx, dx);
    this->fy = fftfreq(ny, dy);
}


std::vector<std::vector<std::complex<double>>> Tensor::peredat_func(double z) const{
    double k = 2.0 * M_PI / lam;
    std::vector<std::vector<std::complex<double>>> H(nx, std::vector<std::complex<double>>(ny));
    for (size_t row = 0; row < nx; ++row){
        double fy_val = fy[row];
        for (size_t col = 0; col < ny; ++col){
            double fx_val = fx[col];
            double arg = 1.0 - (lam * fx_val) * (lam * fx_val) - (lam * fy_val) * (lam * fy_val);
            if (arg < 0.0) arg = 0.0;
            H[row][col] = std::polar(1.0, -k * z * std::sqrt(arg));
        }
    }
    return H;
}


void Tensor::gerchberg_saxton(int iterations, double tolerance){
    const double k = 2.0 * M_PI / lam;
    std::vector<std::vector<std::complex<double>>> field(nx, std::vector<std::complex<double>>(ny));
    for (size_t i = 0; i < nx; ++i)
        for (size_t j = 0; j < ny; ++j)
            field[i][j] = std::polar(std::sqrt(data[0][i][j]), 0.0);
    
    amplitude.assign(nx, std::vector<double>(ny, 0.0));
    phase.assign(nx, std::vector<double>(ny, 0.0));
    errors.clear();
    

    //интерации работы алгоритма
    for(int i = 0; i < iterations; ++i){
        auto prev_amplitude = amplitude;
        

        // прямой проход
        for (size_t i = 1; i < data.size(); ++i) {
            double dz = z_distance[i] - z_distance[i-1];
            auto H = peredat_func(dz);

            auto field_fft = fft2d(field);
            for (size_t r = 0; r < nx; ++r){
                for (size_t c = 0; c < ny; ++c){
                    field_fft[r][c] *= H[r][c];
                }
            }
            field = ifft2d(field_fft);
            for (size_t r = 0; r < nx; ++r){
                for (size_t c = 0; c < ny; ++c) {
                    double amp = std::sqrt(data[i][r][c]);
                    double ph = std::arg(field[r][c]);
                    field[r][c] = std::polar(amp, ph);
                }
            }
        }
        // обратный проход
        for (int i = static_cast<int>(data.size()) - 2; i >= 0; --i) {
            double dz = -z_distance[i+1] + z_distance[i];
            auto H = peredat_func(dz);

            auto field_fft = fft2d(field);

            for (size_t r = 0; r < nx; ++r)
                for (size_t c = 0; c < ny; ++c)
                    field_fft[r][c] *= H[r][c];

            field = ifft2d(field_fft);

            if (i > 0) {
                for (size_t r = 0; r < nx; ++r)
                    for (size_t c = 0; c < ny; ++c) {
                        double amp = std::sqrt(data[i][r][c]);
                        double ph = std::arg(field[r][c]);
                        field[r][c] = std::polar(amp, ph);
                    }
            }
        }

        // Извлекаем амплитуду и фазу на первой плоскости
        for (size_t r = 0; r < nx; ++r) {
            for (size_t c = 0; c < ny; ++c) {
                amplitude[r][c] = std::abs(field[r][c]);
                phase[r][c] = std::arg(field[r][c]);
            }
        }
        // Вычисляем ошибку 
        double error = 0.0;
        for (size_t r = 0; r < nx; ++r)
            for (size_t c = 0; c < ny; ++c)
                error += std::pow(amplitude[r][c] - prev_amplitude[r][c], 2);
        error = std::sqrt(error / (nx * ny));
        errors.push_back(error);

        std::cout << "Итерация" << i+1 << ", ошибка = " << error << std::endl;
        if (error < tolerance) {
            std::cout << "Сходимость достигнута на итерации " << i+1 << std::endl;
            break;
        }
        if (i > 2){
            if ((error > errors[i-1]) && (error > errors[i-2]))
            break;
        }
    }
}