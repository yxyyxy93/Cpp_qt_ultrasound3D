#ifndef UTILS_H
#define UTILS_H

#endif // UTILS_H

#include <complex>
#include <QVector>

QVector<std::complex<double>> analyticSignal(const QVector<double>& signal);

void fft(QVector<std::complex<double>>& signal);

void ifft(QVector<std::complex<double>>& signal);

// surface align utils
template <class T>
void shiftVector_1D(QVector<T>& signal, int k){
    for (int i=0; i<k; ++i){
        signal.append(signal.takeFirst());
    }
}

// 2d analytic-signal
// The first order 2D convolution kernels will be calculated by
double Kernel1(double x, double y, double s);

double Kernel2(double x, double y, double s);

void AnalyticSignal(double cx, double cy,
                    double& Orientation, double& Phase,
                    double& Amplitude, double& ApexAngle,
                    int n, double s_c, double s_f, double dx,
                    QVector<QVector<double>>& img);

QVector<double> min_max_2d(const QVector<QVector<double>> &data);

// To overload the +, -, * operators with the input type QVector<QVector<double>>, point-wise operation
QVector<QVector<double>> operator+(const QVector<QVector<double>>& v1, const QVector<QVector<double>>& v2);
QVector<QVector<double>> operator+(const QVector<QVector<double>>& v1, const double& num);

QVector<QVector<double>> operator-(const QVector<QVector<double>>& v1, const QVector<QVector<double>>& v2);

QVector<QVector<double>> operator*(const QVector<QVector<double>>& v1, const QVector<QVector<double>>& v2);
QVector<QVector<double>> operator*(const QVector<QVector<double>>& v1, const double& num);
