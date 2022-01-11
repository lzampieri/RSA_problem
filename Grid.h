#ifndef THEGRID_H
#define THEGRID_H

#define _HAS_STD_BYTE 0

#include <cstdlib>
#include <cmath>
#include <vector>
#include <tuple>
#include <iostream>
#include <fstream>

#define M_PI 3.14159265

class FourierCoupledGrids;
class GridFiller;

class GridSite {
public:
    static const int Defect = 1;
    static const int Free = 0;
    static const int Atom = -1;
};

template<class T>
class Grid {
protected:
    std::vector<T>* u;

public:
    const int d1, d2;
    const int imax() const { return d1*d2; };

    Grid(int d1, int d2);
    Grid(int d1);
    Grid(const Grid&) = delete; // Copy costruction forbidden
    ~Grid();

    void normalize();
    void gaussian_center_and_normalize();

    void print_data(const char* filename) const;

    // Access
    T& operator[](const int i) const; // i \in [0,d1*d2[
    T& operator()(const int i) const; // i \in [0,d1*d2[
    T& operator()(const std::pair<int,int> xy) const; // x \in [0,d1[, y \in [0,d2[
    T& operator()(const int x, const int y) const; // x \in [0,d1[, y \in [0,d2[

protected:
    // PBC
    const int _pbc_x( int x ) const;
    const int _pbc_y( int y ) const;
    const int _pbc_i( int i ) const;
    const std::pair< int, int > _pbc( const std::pair< int, int > xy ) const;
    const double _abs( std::pair< double, double > xy ) const;

public:
    // Coordinates converting
    const double _i(const int x, const int y) const; // x \in [0,d1[, y \in [0,d2[
    const std::pair< int, int > _xy(const int i) const;  // i \in [0,d1*d2[

public:
    // Distances
    const double d(int i1, int i2) const; // i1,i2 \in [0,d1*d2[
    const double d(int x1, int y1, int x2, int y2) const; // x1,x2 \in [0,d1[; y1,y2 \in [0,d2[
    const double d( std::pair< int, int > xy1, std::pair< int, int > xy2 ) const; // x1,x2 \in [0,d1[; y1,y2 \in [0,d2[

    // Friend classes
    friend class GridFiller;
    friend class FourierCoupledGrids;
};

// Ensure types compilation
template class Grid<double>;
template class Grid<int>;

#endif