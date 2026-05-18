#pragma once
#include <vector>
#include <string>
#include <complex>
#include "fft.h"

class Tensor{
protected:
    std::vector<std::vector<std::vector<double>>> data;
    std::vector<double> z_distance;//расположение плоскостей
    double lam;//длина волны
    std::vector<double> fx;//пространственные частоты(считаем сразу в конструкторе)
    std::vector<double> fy;
    double dx;    //размер пикселя по x
    double dy;    //размер пикселя по y
    size_t nx;    //размер матрицы после препроцессинга(кол-во строк)
    size_t ny;    //размер матрицы после препроцессинга(столбцов)
public:
    std::vector<std::vector<double>> amplitude;
    std::vector<std::vector<double>> phase;
    std::vector<double> errors;



    Tensor(const std::vector<std::string>& path, const std::vector<double> z_distance, double dx, double dy, double lam);
    // Конструктор принимает вектор с путями на все файлы
    // Заполняет тензор интенсивностями
    // Заполняет fx, fy, z


    std::vector<std::vector<std::complex<double>>> peredat_func(double z) const;
    // Считает передаточную функцию на растояние z


    
    void  gerchberg_saxton(int iterations = 100, double tolerance = 1e-4);
    //агоритм gerchberg_saxton {возвращает комплексное поле в 0 плоскости, нужно ли это???}
    // заполняет поля рассчитаными алгоритмом данными amplitude; phase;(в 0 плоскости) errors;

    const std::vector<std::vector<std::vector<double>>>& get_data() const { return data; }
    const std::vector<std::vector<double>>& get_amplitude() const { return amplitude; }
    const std::vector<std::vector<double>>& get_phase() const { return phase; }
    const std::vector<double>& get_errors() const { return errors; }

};