#include "fft.h"
constexpr double PI = 3.14159265358979323846;

void swap(std::vector<std::complex<double>>& vctr, unsigned int i, unsigned int j){
    std::complex<double> buffer;
    buffer = vctr[i];
    vctr[i] = vctr[j];
    vctr[j] = buffer;
}

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

std::vector<std::complex<double>> fftfreq(const int n, const float d = 1){
    std::complex<double> I(0.0, 1.0); 
    std::vector<std::complex<double>> w;
    for (int i = 0; i < n; i+=d){
        w.push_back(std::exp(- std::complex<double>(2 * PI * i/n) * I));
    }
    return w;
}

std::vector<std::complex<double>> ifftfreq(const int n, const float d = 1){
    std::complex<double> I(0.0, 1.0); 
    std::vector<std::complex<double>> w;
    for (int i = 0; i < n; i+=d){
        w.push_back(std::exp(std::complex<double>(2 * PI * i/n) * I));
    }
    return w;
}

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
        std::vector<std::complex<double>> w = fftfreq(m, 1);
        for (int i = 0; i < k; i+=m){
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
        std::vector<std::complex<double>> w = ifftfreq(m, 1);
        for (int i = 0; i < k; i+=m){
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