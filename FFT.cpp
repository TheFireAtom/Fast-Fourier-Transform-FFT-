#include <iostream>
#include <complex>
#include <vector>
#include <cmath>

const double PI = acos(-1);

// Recursive Cooley-Tukey FFT implementation
void fft(std::vector<std::complex<double>>& signal) {
    int N = signal.size();
    if (N <= 1) return;

    // Divide: Split signal into even and odd parts
    std::vector<std::complex<double>> even, odd;
    for (int i = 0; i < N; i += 2) {
        even.push_back(signal[i]);
        odd.push_back(signal[i + 1]);
    }

    // Recursive FFT for even and odd parts
    fft(even);
    fft(odd);

    // Conquer: Combine results
    for (int k = 0; k < N / 2; ++k) {
        std::complex<double> t = std::polar(1.0, -2 * PI * k / N) * odd[k];
        signal[k] = even[k] + t;
        signal[k + N / 2] = even[k] - t;
    }
}

// Generate input signal using specified function
void generateSignal(std::vector<double>& signal, double fs) {
    const double ts = 1.0 / fs;
    for (size_t n = 0; n < signal.size(); ++n) {
        signal[n] = 1.5 * std::cos(2 * PI * 1500 * n * ts + (PI / 4)) +
            0.25 * std::cos(2 * PI * 1600 * n * ts + (PI / 2)) +
            1 * std::cos(2 * PI * 4500 * n * ts + (PI / 4));
    }
}

int main() {
    const int N = 4096;
    const double fs = 18000.0;

    std::vector<double> signal(N);
    generateSignal(signal, fs);

    std::vector<std::complex<double>> complex_signal(N);
    for (int i = 0; i < N; ++i) {
        complex_signal[i] = signal[i];
    }

    fft(complex_signal);

    // Print FFT result...
    for (int i = 0; i < N; ++i) {
        std::cout << "FFT[" << i << "] = " << complex_signal[i] << std::endl;
    }

    return 0;
}
