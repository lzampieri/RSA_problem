#ifndef THEGRID_H
#define THEGRID_H

#include <cstdlib>
#include <cmath>
#include <vector>
#include <ctime>
#include <tuple>
#include <iostream>
#include <fstream>
#include "fftw3.h"
#include "CorrFuncCalcolator.h"

using namespace std;

class Grid {
protected:
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

    // Access
    double operator()(const int i) const; // i \in [0,d1*d2[
    double operator()(const pair<int,int> xy, double fourier = false) const; // x \in [0,d1[, y \in [0,d2[
    double operator()(const int x, const int y, double fourier = false) const; // x \in [0,d1[, y \in [0,d2[

protected:
    // PBC
    const int _pbc_x( const int x ) const;
    const int _pbc_y( const int y ) const;
    const int _pbc_i( const int i ) const;
    const pair< int, int > _pbc( const pair< int, int > xy ) const;

    // Coordinates converting
    const double _i(const int x, const int y) const; // x \in [0,d1[, y \in [0,d2[
    const pair< int, int > _xy(const int i) const;  // i \in [0,d1*d2[

public:
    // Distances
    const double d(int i1, int i2) const; // i1,i2 \in [0,d1*d2[
    const double d(int x1, int y1, int x2, int y2) const; // x1,x2 \in [0,d1[; y1,y2 \in [0,d2[
    const double d( pair< int, int > xy1, pair< int, int > xy2 ) const; // x1,x2 \in [0,d1[; y1,y2 \in [0,d2[

    // Friend classes
    friend class CorrFuncCalcolator;
};

#endif