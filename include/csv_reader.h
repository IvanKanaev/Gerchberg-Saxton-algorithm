#pragma once

#include <string>
#include <vector>
#include "types.h"


Matrix_2D read_csv_to_matrix(const std::string& filepath, char delimiter = ';');
// Читает CSV-файл с интенсивностями данными.
// @param filepath   путь к CSV-файлу
// @param delimiter  символ-разделитель (по умолчанию ';')
// @return           матрица типа Matri2D (vector<vector<double>>)
