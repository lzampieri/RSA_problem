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
    int d1d2;
    int imax() const { return d1*d2; };

    GridProps(int d1, int d2);
    GridProps(int d);
    GridProps(const GridProps& G);

protected:
    // PBC
    inline int _pbc_x( int x ) const { return x < 0 ? x + d1 : x % d1; }
    inline int _pbc_y( int y ) const { return y < 0 ? y + d2 : y % d2; }
    inline int _pbc_i( int i ) const { return i < 0 ? i + d1d2 : i % d1d2; }

public:
    // Coordinates converting
    inline double _i(const int x, const int y) const { return _pbc_x(x) * d2 + _pbc_y(y); } // x \in [0,d1[, y \in [0,d2[
    inline GridSite _xy(const int i) const;  // i \in [0,d1*d2[

public:
    // Distances
    inline double d(int i1, int i2) const; // i1,i2 \in [0,d1*d2[
    double d(int x1, int y1, int x2, int y2) const; // x1,x2 \in [0,d1[; y1,y2 \in [0,d2[
    inline double d(const GridSite& xy1, const GridSite& xy2 ) const; // x1,x2 \in [0,d1[; y1,y2 \in [0,d2[

    // Inline functions declared lower in this file
    // 'cause they need the complete implementation of GridSite
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
    inline int I() const { return _i(X, Y); }

    // Operators
    inline void setX(int newX) { _X = _pbc_x(newX); }
    inline void setY(int newY) { _Y = _pbc_y(newY); }
    GridSite& operator= (const GridSite& gs);
    GridSite& operator+=(const GridSite& xy2);
    inline friend GridSite operator+(const GridSite& xy, const GridSite& xy2) { return GridSite(xy.X + xy2.X, xy.Y + xy2.Y, xy); }
    GridSite& operator-=(const GridSite& xy2);
    inline friend GridSite operator-(const GridSite& xy, const GridSite& xy2) { return GridSite(xy.X - xy2.X, xy.Y - xy2.Y, xy); }
    inline friend bool operator== (const GridSite& xy1, const GridSite& xy2) { return xy1.X == xy2.X && xy1.Y == xy2.Y; }
    inline friend bool operator!= (const GridSite& xy1, const GridSite& xy2) { return !(xy1 == xy2); }
    inline double abs() const { return sqrt((double)X * X + (double)Y * Y); }

    // Distances
    using GridProps::d;
    inline double d(const GridSite& xy2) const { return d(*this, xy2); }

    // Utils
    inline bool last_row() const { return _Y == d2 - 1; }
    
    static const short Defect = -1;
    static const short Free = 0;
    static const short Atom = 1;
};

inline GridSite GridProps::_xy(const int i) const { return GridSite(_pbc_i(i) / d2, _pbc_i(i) % d2, *this); }
inline double GridProps::d(int i1, int i2) const { return d(_xy(i1), _xy(i2)); }; // i1,i2 \in [0,d1*d2[
inline double GridProps::d(const GridSite& xy1, const GridSite& xy2 ) const { return d(xy1.X, xy1.Y, xy2.X, xy2.Y); }; // x1,x2 \in [0,d1[; y1,y2 \in [0,d2[

template<class T>
class Grid : public GridProps {
protected:
    std::vector<T>* u;

public:
    Grid(int d1, int d2);
    Grid(int d);
    Grid(const Grid&) = delete; // Copy costruction forbidden
    ~Grid();

    void print_data(const std::string filename) const;

    // Access
    inline T& operator[](const int i) const { return u->at(_pbc_i(i)); }; // i \in [0,d1*d2[
    inline T& operator[](const GridSite xy) const { return u->at(xy.I()); }; // x \in [0,d1[, y \in [0,d2[
    inline T& operator()(const int x, const int y) const { return operator[](GridSite(x, y, *this)); }; // x \in [0,d1[, y \in [0,d2[
    inline T& at(const int i) const { return operator[](i); };
    inline T& at(const GridSite xy) const { return operator[](xy); };
    inline T& at(const int x, const int y) const { return operator()(x,y); };

public:
    // Friend classes
    friend class GridFiller;
    friend class FourierGrids;
};

// Ensure types compilation
template class Grid<double>;
template class Grid<int>;

#endif