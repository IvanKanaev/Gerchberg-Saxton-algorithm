#include "tensor.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>


int main() {
    const int min_z = -7;
    const int max_z = 13;
    const int num_planes = 21;
    double weight = 0.02;
    double gain = 1.0;

    const std::string folder = "full_real_data/";
    
    std::vector<std::string> paths(num_planes);
    std::vector<double> z_distance(num_planes);

    for (int i = min_z; i <= max_z; ++i) {
        int idx = i - min_z; 
        paths[idx] = folder + "full.z=" + std::to_string(i) + "mm.csv";
        z_distance[idx] = i * 1e-3; 
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

    t.save_amplitude("amplitude.csv");
    t.save_phase("phase.csv");
    return 0;
}