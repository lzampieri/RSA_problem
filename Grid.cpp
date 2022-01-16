#include "Grid.h"

using namespace std;

// GridProps
// Constructors
GridProps::GridProps(int d1, int d2) : d1(d1), d2(d2) {}
GridProps::GridProps(int d) : d1(d), d2(d) {}
GridProps::GridProps(const GridProps& G) : d1(G.d1), d2(G.d2) {}

// PBC
int GridProps::_pbc_x( int x ) const {
    while( x < 0 ) x += d1;
    return x % d1;
}

int GridProps::_pbc_y( int y ) const {
    while( y < 0 ) y += d2;
    return y % d2;
}

int GridProps::_pbc_i( int i ) const {
    while( i < 0 ) i += d1*d2;
    return i % (d1*d2);
}

// Coordinates converting

double GridProps::_i(const int x, const int y) const {
    return _pbc_x(x) * d2 + _pbc_y(y);
}

GridSite GridProps::_xy(const int i) const {
    return GridSite( _pbc_i(i) / d2, _pbc_i(i) % d2, *this );   
}

// Distances

double GridProps::d(int i1, int i2) const{
    return d( _xy(i1), _xy(i2) );
}

double GridProps::d(int x1, int y1, int x2, int y2) const {
    double deltaX = _pbc_x( x1 - x2 );
    deltaX = min( deltaX, d1 - deltaX );
    double deltaY = _pbc_y( y1 - y2 );
    deltaY = min( deltaY, d2 - deltaY );
    return sqrt( deltaX * deltaX + deltaY * deltaY );
}

double GridProps::d(const GridSite& xy1, const GridSite& xy2 ) const {
    return d( xy1.X, xy1.Y, xy2.X, xy2.Y );
}

// GridSite
// Constructors

GridSite::GridSite(int X, int Y, const GridProps& gp) : GridProps(gp) { setX(X); setY(Y); }
GridSite::GridSite(int I, const GridProps& gp) : GridSite( gp._xy(I) ) {}
GridSite::GridSite(const GridSite& G) : GridProps( G ) { setX(G.X); setY(G.Y); }

// Access to I
int GridSite::I() const {
    return _i( X, Y );
}

// Operators

void GridSite::setX(int newX) {
    _X = _pbc_x( newX );
}

void GridSite::setY(int newY) {
    _Y = _pbc_y( newY );
}

GridSite& GridSite::operator= (const GridSite& gs) {
    d1 = gs.d1;
    d2 = gs.d2;
    setX( gs.X );
    setY( gs.Y );
    return *this;
}

GridSite& GridSite::operator+= (const GridSite& xy2) {
    setX( X + xy2.X );
    setY( Y + xy2.Y );
    return *this;
}

GridSite operator+ (const GridSite& xy, const GridSite& xy2) {
    return GridSite( xy.X + xy2.X, xy.Y + xy2.Y, xy );
}

GridSite& GridSite::operator-=(const GridSite& xy2) {
    setX( X - xy2.X );
    setY( Y - xy2.Y );
    return *this;
}

GridSite operator-(const GridSite& xy, const GridSite& xy2) {
    return GridSite( xy.X - xy2.X, xy.Y - xy2.Y, xy );
}

bool operator== (const GridSite& xy1, const GridSite& xy2) {
    return xy1.X == xy2.X && xy1.Y == xy2.Y;
}

bool operator!= (const GridSite& xy1, const GridSite& xy2) {
    return !( xy1 == xy2 );
}

double GridSite::abs() const {
    return sqrt( (double)X * X + (double)Y * Y );
}

// Distances

double GridSite::d(const GridSite& xy2) {
    return d(*this, xy2);
}

// Grid
// Constructors

template<class T>
Grid<T>::Grid(int d) : Grid<T>(d, d) {};

template<class T>
Grid<T>::Grid(int d1, int d2) : GridProps(d1,d2) {
    u = new vector<T>(d1*d2, 0);
}

template<class T>
Grid<T>::~Grid() {
    u->clear();
    delete u;
}

template<class T>
void Grid<T>::normalize() {
    T maximum = u->at(0);
    for(int i=1; i<d1*d2; i++) {
        maximum = max( maximum, u->at(i) );
    }
    for(int i=0; i<d1*d2; i++) {
        u->at(i) /= maximum;
    }
}

template<class T>
void Grid<T>::gaussian_center_and_normalize() {
    double sum = 0;
    double sum2 = 0;
    for(int i=0; i<d1*d2; i++) {
        sum += u->at(i);
        sum2+= u->at(i)*u->at(i);
    }
    double avg = sum / (d1*d2);
    double sigma = sqrt( ( sum2 - pow( sum, 2 ) / (d1*d2) ) / ( d1*d2 - 1 ) );

    for(int i=0; i<d1*d2; i++) {
        u->at(i) = ( u->at(i) - avg ) / sigma;
    }
}

template<class T>
void Grid<T>::print_data(const char* filename) const {
    ofstream out(filename);
    for(int i=0; i<d1; i++)
        for(int j=0; j<d2; j++) {
            out<<i<<'\t'<<j<<'\t'<<operator()(i,j)<<'\n';
        }
    out.close();
}

// Access

template<class T>
T& Grid<T>::operator[](const int i) const {
    return operator()(i);
}

template<class T>
T& Grid<T>::operator()(const int i) const {
    return u->at( _pbc_i(i) );
}

template<class T>
T& Grid<T>::operator()(const GridSite xy) const {
    return u->at( xy.I() );
}

template<class T>
T& Grid<T>::operator()(const int x, const int y) const {
    return operator()( GridSite( x, y, *this ) );
}