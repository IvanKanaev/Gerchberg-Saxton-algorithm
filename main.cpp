#include "tensor.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>


int main()
{
    std::vector<std::string> paths(18);
    for (int i = 0; i < 18; ++i) {
        paths[i] = "test_data/" + std::to_string(i + 1) + ".csv";
    }
    std::vector<double> z_distance(paths.size());
    for (size_t i = 0; i < paths.size(); ++i){
        z_distance[i] = i * 0.5e-3;
    }   
    double dx = 1.85e-6;
    double dy = 1.85e-6;
    double lam = 660e-9;
    const int max_iter = 100;
    const double tolerance = 1e-6;  

    Tensor t(paths, z_distance, dx, dy, lam);
    t.gerchberg_saxton(max_iter, tolerance); 
    const auto& errors = t.get_errors();
    std::cout << "Итоговая ошибка: " << errors.back() << " после " << errors.size() << " итерации\n";
    for (size_t i = 0; i < errors.size(); ++i){
        std::cout << "итерация  " << i+1 << ": " << errors[i] << "\n";
    }

    return 0;
}