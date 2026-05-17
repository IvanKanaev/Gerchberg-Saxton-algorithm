#include "tensor.h"
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