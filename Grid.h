#include <cstdlib>
#include <cmath>
#include <vector>
#include <ctime>
#include <tuple>
#include <iostream>
#include <fstream>
#include "fftw3.h"


using namespace std;

class Grid {
    const int d1, d2;
    vector<double>* h;
    vector<double>* f;
    fftw_plan direct_plan, reverse_plan;

public:
    Grid(int d1, int d2);
    Grid(int d1);

    void populate_iid();
    void populate_coscos(double fact1, double fact2);
    void populate_square(int side);

    void fourier_transform(double reverse = false);

    void multiply_fft(double gamma);

    void print_data(const char* filename, double fourier = false) const;

    double operator()(const int i) const;
    double operator()(const int i, const int j, double fourier = false) const;

    const double d(int i1, int i2) const;
    const double d(int x1, int y1, int x2, int y2) const;

    vector< tuple<double,double, double> >* corr_func(double resol) const;
    void print_corr(const char* filename) const;
};