#include "utils.h"

#include <cmath>
#include <complex>
#include <QVector>
#include <fftw3.h>

// analytic-signal calculation
QVector<std::complex<double>> analyticSignal(const QVector<double>& signal)
{
    int N = signal.size();

    // Compute the FFT of the signal
    QVector<std::complex<double>> signalFFT(N);
    for (int i = 0; i < N; i++) {
        signalFFT[i] = std::complex<double>(signal[i], 0.0);
    }
    fft(signalFFT); // Replace with your own FFT function

    // Compute the frequency domain representation of the Hilbert transform
    QVector<std::complex<double>> hilbertFFT(N);
    hilbertFFT[0] = std::complex<double>(0.0, 0.0);
    for (int i = 1; i < N/2; i++) {
        hilbertFFT[i] = std::complex<double>(0.0, -1.0/(M_PI*i));
        hilbertFFT[N-i] = std::conj(hilbertFFT[i]);
    }
    hilbertFFT[N/2] = std::complex<double>(0.0, 0.0);

    // Multiply the FFT of the signal by the frequency domain representation of the Hilbert transform
    QVector<std::complex<double>> analyticFFT(N);
    for (int i = 0; i < N/2; i++) {
        analyticFFT[i] = signalFFT[i] * std::complex<double>(2.0, 0.0);
    }
    for (int i = N/2; i < N; i++) {
        analyticFFT[i] = std::complex<double>(0.0, 0.0);
    }

    // Compute the inverse FFT of the analytic signal
    ifft(analyticFFT); // Replace with your own IFFT function

    return analyticFFT;
}

void fft(QVector<std::complex<double>>& signal)
{
    int N = signal.size();

    // Create a FFTW plan for forward FFT
    fftw_plan plan = fftw_plan_dft_1d(N, reinterpret_cast<fftw_complex*>(&signal[0]),
            reinterpret_cast<fftw_complex*>(&signal[0]),
            FFTW_FORWARD, FFTW_ESTIMATE);

    // Execute the FFTW plan
    fftw_execute(plan);

    // Destroy the FFTW plan
    fftw_destroy_plan(plan);
}

void ifft(QVector<std::complex<double>>& signal)
{
    int N = signal.size();

    // Create a FFTW plan for backward FFT
    fftw_plan plan = fftw_plan_dft_1d(N, reinterpret_cast<fftw_complex*>(&signal[0]),
            reinterpret_cast<fftw_complex*>(&signal[0]),
            FFTW_BACKWARD, FFTW_ESTIMATE);

    // Execute the FFTW plan
    fftw_execute(plan);

    // Normalize the output of IFFT
    for (int i = 0; i < N; i++) {
        signal[i] /= N;
    }

    // Destroy the FFTW plan
    fftw_destroy_plan(plan);
}

// 2d analytic-signal
// The first order 2D convolution kernels will be calculated by
double Kernel1(double x, double y, double s)
{
    double ss = pow(s, 2);
    double kk = pow(x, 2) + pow(y, 2);
    return 1 / (2 * M_PI * pow(ss + kk, 1.5));
}
//the second order 2D convolution kernels in spatial domain(x, y) with scale space parameter s will be determined by

double Kernel2(double x, double y, double s)
{
    double ss = pow(s, 2);
    double kk = pow(x, 2) + pow(y, 2);
    double d = pow(kk, 2) * pow(ss + kk, 1.5) * 2 * M_PI;
    return (d == 0) ? 0 : (s * (2 * ss + 3 * kk) - 2 * pow(ss + kk, 1.5)) / d;
}

void AnalyticSignal(double cx, double cy,
                    double& Orientation, double& Phase,
                    double& Amplitude, double& ApexAngle,
                    int n, double s_c, double s_f, double dx,
                    QVector<QVector<double>>& img)
{
    double f_p = 0, f_x = 0, f_y = 0, f_xx = 0, f_xy = 0, f_yy = 0;
    //2D convolution
    for (int x = -n;x <= n;x += dx)
        for (int y = -n;y <= n;y += dx)
        {
            double t = img[x + cx][y + cy] * pow(dx, 2);
            double pf = t * Kernel1(x, y, s_f);
            double pc = t * Kernel1(x, y, s_c);
            double k = t * (Kernel2(x, y, s_f) - Kernel2(x, y, s_c));
            //signal in Poisson scale space
            f_p += s_f * pf - s_c * pc;
            //first order Hilbert transform
            f_x += x * (pf - pc);
            f_y += y * (pf - pc);
            //second order Hilbert transform
            f_xx += x * x * k;
            f_yy += y * y * k;
            f_xy += x * y * k;
        }
    double f_pm = 0.5 * (f_xx - f_yy);
    double f_s = 0.5 * f_p;
    double f_plus = f_xy;
    double e = sqrt(pow(f_pm, 2) + pow(f_plus, 2)) / fabs(f_s);
    double q = (pow(f_x, 2) + pow(f_y, 2)) * 2 / (1 + e);
    Phase = atan2(sqrt(q), f_p);
    Orientation = ((int)Phase == 0)
            ? 0.5 * atan2(f_plus, f_pm) + M_PI / 2
            : atan2(f_y, f_x);
    Amplitude = 0.5 * sqrt(pow(f_p, 2) + q);
    ApexAngle = atan2(sqrt(pow(f_s, 2) - pow(f_plus, 2) -
                           pow(f_pm, 2)), sqrt(pow(f_plus, 2) + pow(f_pm, 2)));
}

QVector<double> min_max_2d(const QVector<QVector<double>> &data){
    // initialize the minimum and maximum values to the first element in the vector
    double minVal = data[0][0];
    double maxVal = data[0][0];

    // iterate over all the elements in the vector and update the minimum and maximum values
    for (int i = 0; i < data.size(); ++i) {
        for (int j = 0; j < data[i].size(); ++j) {
            if (data[i][j] < minVal) {
                minVal = data[i][j];
            }
            if (data[i][j] > maxVal) {
                maxVal = data[i][j];
            }
        }
    }
    QVector<double> results;
    results.push_back(minVal);
    results.push_back(maxVal);

    return results;
}


// To overload the +, -, * operators with the input type QVector<QVector<double>>, point-wise operation
QVector<QVector<double>> operator+(const QVector<QVector<double>>& v1, const QVector<QVector<double>>& v2) {
    int numRows = v1.size();
    int numCols = v1[0].size();

    QVector<QVector<double>> result(numRows, QVector<double>(numCols));

    for (int i = 0; i < numRows; i++) {
        for (int j = 0; j < numCols; j++) {
            result[i][j] = v1[i][j] + v2[i][j];
        }
    }

    return result;
}
QVector<QVector<double>> operator+(const QVector<QVector<double>>& v1, const double& num){
    int numRows = v1.size();
    int numCols = v1[0].size();

    QVector<QVector<double>> result(numRows, QVector<double>(numCols));

    for (int i = 0; i < numRows; i++) {
        for (int j = 0; j < numCols; j++) {
            result[i][j] = v1[i][j] + num;
        }
    }

    return result;
}

QVector<QVector<double>> operator-(const QVector<QVector<double>>& v1, const QVector<QVector<double>>& v2) {
    int numRows = v1.size();
    int numCols = v1[0].size();

    QVector<QVector<double>> result(numRows, QVector<double>(numCols));

    for (int i = 0; i < numRows; i++) {
        for (int j = 0; j < numCols; j++) {
            result[i][j] = v1[i][j] - v2[i][j];
        }
    }

    return result;
}

QVector<QVector<double>> operator*(const QVector<QVector<double>>& v1, const QVector<QVector<double>>& v2) {
    int numRows = v1.size();
    int numCols = v1[0].size();

    QVector<QVector<double>> result(numRows, QVector<double>(numCols));

    for (int i = 0; i < numRows; i++) {
        for (int j = 0; j < numCols; j++) {
            result[i][j] = v1[i][j] * v2[i][j];
        }
    }

    return result;
}
QVector<QVector<double>> operator*(const QVector<QVector<double>>& v1, const double& num) {
    int numRows = v1.size();
    int numCols = v1[0].size();

    QVector<QVector<double>> result(numRows, QVector<double>(numCols));

    for (int i = 0; i < numRows; i++) {
        for (int j = 0; j < numCols; j++) {
            result[i][j] = v1[i][j] * num;
        }
    }

    return result;
}
