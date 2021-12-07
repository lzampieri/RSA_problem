#include "Grid.h"

using namespace std;

template<class T>
Grid<T>::Grid(int d1) : Grid(d1, d1) {};

template<class T>
Grid<T>::Grid(int d1, int d2) : d1(d1), d2(d2) {
    u = new vector<T>(d1*d2, 0);
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
        for(int j=0; j<d2; j++)
            out<<i<<'\t'<<j<<'\t'<<operator()(i,j)<<'\n';
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
T& Grid<T>::operator()(const pair<int,int> xy) const {
    return operator()( xy.first, xy.second );
}

template<class T>
T& Grid<T>::operator()(const int x, const int y) const {
    return u->at( _i(x,y) );
}

// PBC

template<class T>
const int Grid<T>::_pbc_x( int x ) const {
    while( x < 0 ) x += d1;
    return x % d1;
}

template<class T>
const int Grid<T>::_pbc_y( int y ) const {
    while( y < 0 ) y += d1;
    return y % d2;
}

template<class T>
const int Grid<T>::_pbc_i( int i ) const {
    while( i < 0 ) i += d1*d2;
    return i % (d1*d2);
}

template<class T>
const pair< int, int > Grid<T>::_pbc( const pair< int, int > xy ) const {
    return make_pair( _pbc_x( xy.first ), _pbc_y( xy.second ) );
}

template<class T>
const double Grid<T>::_abs( pair< double, double > xy ) const {
    return sqrt( xy.first * xy.first + xy.second * xy.second );
}

// Coordinates converting
template<class T>
const double Grid<T>::_i(const int x, const int y) const {
    return _pbc_x(x) * d2 + _pbc_y(y);
}

template<class T>
const pair< int, int > Grid<T>::_xy(const int i) const {
    return make_pair< int, int >( _pbc_i(i) / d2, _pbc_i(i) % d2 );   
}

// Distances
template<class T>
const double Grid<T>::d(int i1, int i2) const{
    return d( _xy(i1), _xy(i2) );
}

template<class T>
const double Grid<T>::d(int x1, int y1, int x2, int y2) const{
    int deltaX = _pbc_x( x1 - x2 );
    deltaX = min( deltaX, d1 - deltaX );
    int deltaY = _pbc_y( y1 - y2 );
    deltaY = min( deltaY, d2 - deltaY );
    return sqrt( deltaX * deltaX + deltaY * deltaY );
}

template<class T>
const double Grid<T>::d( pair< int, int > xy1, pair< int, int > xy2 ) const {
    return d( xy1.first, xy1.second, xy2.first, xy2.second );
}