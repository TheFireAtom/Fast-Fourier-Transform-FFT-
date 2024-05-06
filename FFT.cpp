// ����������� ����������� ���������
#include <iostream>
#include <complex>
#include <vector>
#include <cmath>

const double PI = acos(-1); // ��������� �� ������� � ������� ����������� �� -1

// ��������� ��������� ���� � ����� � ������� ��������
void fft(std::vector<std::complex<double>>& signal) {
    int N = signal.size();
    if (N <= 1) return;

    // ���������� ��������� ������� �� ������ � �������� �����
    std::vector<std::complex<double>> even, odd;
    for (int i = 0; i < N; i += 2) {
        even.push_back(signal[i]);
        odd.push_back(signal[i + 1]);
    }

    // �������� ��� ��� ������ � �������� �����
    fft(even);
    fft(odd);

    // ���������� �����������
    for (int k = 0; k < N / 2; ++k) {
        // ������� polar ������� ����������� ����� �� ������ �������� ���������, � �������� ���������� ����������� �������� ��������� � ����
        std::complex<double> t = std::polar(1.0, -2 * PI * k / N) * odd[k]; 
        // ��������� �������� "�������", ������� ��������� ��� ������ �������������� ������ ���������
        signal[k] = even[k] + t;
        signal[k + N / 2] = even[k] - t;
    }
}

// �������� ��������� ������� �� ������ ����������� �������
void generateSignal(std::vector<double>& signal, double fs) {
    const double ts = 1.0 / fs;
    for (size_t n = 0; n < signal.size(); ++n) {
        signal[n] = 1.5 * std::cos(2 * PI * 1500 * n * ts + (PI / 4)) +
            0.25 * std::cos(2 * PI * 1600 * n * ts + (PI / 2)) +
            1 * std::cos(2 * PI * 4500 * n * ts + (PI / 4));
    }
}


// �������� ������� � ��������� ��������� ������� 
int main() {
    // �������� N = 4096 ���� ������� ����� ������������������ �������� ������ ���������
    // ��� ��������� ��������� ����������, ����� ���������� �������� ��������������� ������� ������
    const int N = 4096; // ���������� �������� 
    const double fs = 18000.0; // ������� �������������
     
    // �������� ������� (������������� �������) � �������� � N ������
    std::vector<double> signal(N);
    // ����� ������� �������� �������, � �������� ���������� ������ ������ signal � ������� ������������� fs
    generateSignal(signal, fs);

    // ���������������� �������� ������� signal � ������ � ������������ ���������� complex_signal
    std::vector<std::complex<double>> complex_signal(N);
    for (int i = 0; i < N; ++i) {
        complex_signal[i] = signal[i];
    }

    // ����� ������� ��� ��� �������� �������� � ������� complex_signal
    fft(complex_signal);

    // ����� ���������� �����������
    for (int i = 0; i < N; ++i) {
        std::cout << "FFT[" << i << "] = " << complex_signal[i] << std::endl;
    }

    return 0;
}
