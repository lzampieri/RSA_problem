#ifndef THEGRID_H
#define THEGRID_H

#define _HAS_STD_BYTE 0

#include <cstdlib>
#include <cmath>
#include <vector>
#include <tuple>
#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <thread>

#ifndef M_PI
#define M_PI 3.14159265
#endif

class FourierGrids;
class GridFiller;
class GridSite;

class GridProps {
public:
    int d1, d2;
    int imax() const { return d1*d2; };

    GridProps(int d1, int d2);
    GridProps(int d);
    GridProps(const GridProps& G);

protected:
    // PBC
    int _pbc_x( int x ) const;
    int _pbc_y( int y ) const;
    int _pbc_i( int i ) const;

public:
    // Coordinates converting
    double _i(const int x, const int y) const; // x \in [0,d1[, y \in [0,d2[
    GridSite _xy(const int i) const;  // i \in [0,d1*d2[

public:
    // Distances
    double d(int i1, int i2) const; // i1,i2 \in [0,d1*d2[
    double d(int x1, int y1, int x2, int y2) const; // x1,x2 \in [0,d1[; y1,y2 \in [0,d2[
    double d(const GridSite& xy1, const GridSite& xy2 ) const; // x1,x2 \in [0,d1[; y1,y2 \in [0,d2[
};

class GridSite : public GridProps {
private:
    int _X, _Y;

public:
    GridSite(int X, int Y, const GridProps& gp);
    GridSite(int I, const GridProps& gp);
    GridSite(const GridSite& G);

    // Read-only access to X and Y
    const int& X = _X;
    const int& Y = _Y;

    // Access to I
    int I() const;

    // Operators
    void setX(int newX);
    void setY(int newY);
    GridSite& operator= (const GridSite& gs);
    GridSite& operator+=(const GridSite& xy2);
    friend GridSite operator+(const GridSite& xy, const GridSite& xy2);
    GridSite& operator-=(const GridSite& xy2);
    friend GridSite operator-(const GridSite& xy, const GridSite& xy2);
    friend bool operator== (const GridSite& xy1, const GridSite& xy2);
    friend bool operator!= (const GridSite& xy1, const GridSite& xy2);
    double abs() const;

    // Distances
    using GridProps::d;
    double d(const GridSite& xy2);

    // Utils
    bool last_row() const;
    
    static const short Defect = -1;
    static const short Free = 0;
    static const short Atom = 1;
};

template<class T>
class Grid : public GridProps {
protected:
    std::vector<T>* u;

public:
    Grid(int d1, int d2);
    Grid(int d);
    Grid(const Grid&) = delete; // Copy costruction forbidden
    ~Grid();

    void normalize();
    void gaussian_center_and_normalize();

    void print_data(const std::string filename) const;

    // Access
    T& operator[](const int i) const; // i \in [0,d1*d2[
    T& operator[](const GridSite xy) const; // x \in [0,d1[, y \in [0,d2[
    T& operator()(const int x, const int y) const; // x \in [0,d1[, y \in [0,d2[

public:
    // Friend classes
    friend class GridFiller;
    friend class FourierGrids;
};

// Ensure types compilation
template class Grid<double>;
template class Grid<int>;

#endif