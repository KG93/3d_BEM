#ifndef BEMOPERATORSUBROUTINES_H
#define BEMOPERATORSUBROUTINES_H
#include <cmath>
#include <complex>

#include "global.h"

//  "A collocation BEM for 3D acoustic problems based on a non-singular
//  Burton-Miller formulation with linear continuous elements" by Haijun Wu, Wenjing Ye and Weikang Jiang

//class bemOperatorSubroutines
//{
//public:
//    bemOperatorSubroutines();


//};

static inline double l0(int m , double angle)
{
    switch (m)
    {
        case -1:
            return std::sin(angle);

        case 0:
            return angle;

        case 1:
            {
                double sinAngle = std::sin(angle);
                return std::log((1+sinAngle) / (1-sinAngle)) / 2.0;
            }
        default:
            if(m < -1)
            {
                return 0;
            }
            else
            {
                double secAngleM = std::pow(1.0 / std::cos(angle), m-1);
                return std::sin(angle) * secAngleM / (m -1) + (m-2.0) / (m-1.0) * l0(m-2, angle);
            }
    }
}

static inline double l1(int m , double angle)
{
    switch (m)
    {
        case 1:
            return -std::log( std::cos(angle) );

        default:
            return (std::pow(std::cos(angle), 1-m) - 1)/ (m - 1.0);
    }
}

static inline std::complex<double> integ_exp_ikr(double angle, double h, std::complex<double> iWavenumber)
{
    std::complex<double> iWavenumberH =  iWavenumber * h;;
    std::complex<double> tmp = 1;

    int m = 0;
    std::complex<double> result = l0(m, angle);
    while(std::abs(tmp) > global::tiny)
    {
        m++;
        tmp *= iWavenumberH / (double)m;
        result += tmp * l0(m, angle);
    }
    return result;
}

static inline std::complex<double> integ_exp_ikr_R(double angle, double h, std::complex<double> iWavenumber)
{
    std::complex<double> iWavenumberH =  iWavenumber * h;;
    std::complex<double> tmp = 1;

    int m = 0;
    std::complex<double> result = l0(m-1, angle);
    while(std::abs(tmp) > global::tiny)
    {
        m++;
        tmp *= iWavenumberH / (double)m;
        result += tmp * l0(m-1, angle);
    }
    return result/h;
}

static inline std::complex<double> integ_sin_exp_ikr_R(double angle, double h, std::complex<double> iWavenumber)
{
    std::complex<double> iWavenumberH =  iWavenumber * h;;
    std::complex<double> tmp = 1;

    int m = 0;
    std::complex<double> result = l1(m, angle);
    while(std::abs(tmp) > global::tiny)
    {
        m++;
        tmp *= iWavenumberH / (double)m;
        result += tmp * l1(m, angle);
    }
    return result;
}

static inline std::complex<double> integ_cos_exp_ikr_R(double angle, double h, std::complex<double> iWavenumber)
{
    return integ_exp_ikr_R(angle, h, iWavenumber) * h;
}

static inline std::complex<double> integ_sin_R_exp_ikr(double angle, double h, std::complex<double> iWavenumber)
{
    std::complex<double> iWavenumberH =  iWavenumber * h;;
    std::complex<double> tmp = 1;

    int m = 0;
    std::complex<double> result = l1(m+1, angle);
    while(std::abs(tmp) > global::tiny)
    {
        m++;
        tmp *= iWavenumberH / (double)m;
        result += tmp * l1(m+1, angle);
    }
    return result * h;
}

static inline std::complex<double> integ_cos_R_exp_ikr(double angle, double h, std::complex<double> iWavenumber)
{
    return integ_exp_ikr(angle, h, iWavenumber) * h;
}

static inline std::complex<double> integ_ikr_sin(double angle, double h, std::complex<double> iWavenumber)
{
    std::complex<double> iWavenumberH =  iWavenumber * h;;
    std::complex<double> tmp = 1;

    int m = 0;
    std::complex<double> result = 0;
    while(std::abs(tmp) > global::tiny)
    {
        m++;
        tmp *= iWavenumberH / (double)m;
        result += tmp * l1(m, angle) / (double)m;
    }
    return result;
}

static inline std::complex<double> integ_ikr_cos(double angle, double h, std::complex<double> iWavenumber)
{
    std::complex<double> iWavenumberH =  iWavenumber * h;;
    std::complex<double> tmp = 1;

    int m = 0;
    std::complex<double> result = 0;
    while(std::abs(tmp) > global::tiny)
    {
        m++;
        tmp *= iWavenumberH / (double)m;
        result += tmp * l0(m-1, angle) / (double)m;
    }
    return result;
}

static inline double sinLnR(double angle, double h)
{
    double cosAngle = std::cos(angle);
    return cosAngle * std::log(cosAngle) - (1 + std::log(h)) * (cosAngle - 1);
}

static inline double cosLnR(double angle, double h)
{
    double sinAngle = std::sin(angle);
    double cosAngle = std::cos(angle);
    return sinAngle * (std::log(h) - 1) + (std::log((1 + sinAngle) / (1 - sinAngle)) + 2.0 * sinAngle * std::log(cosAngle)) / 2.0;
}
#endif // BEMOPERATORSUBROUTINES_H
