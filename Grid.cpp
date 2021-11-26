#include "Grid.h"

Grid::Grid(int d1) : Grid(d1, d1) {};

Grid::Grid(int d1, int d2) : d1(d1), d2(d2) {
    h = new vector<double>(d1*d2);
    f = new vector<double>(d1*d2);
    direct_plan = fftw_plan_r2r_2d(d1,d2,h->data(),f->data(), FFTW_REDFT00, FFTW_REDFT00, FFTW_ESTIMATE);
    reverse_plan= fftw_plan_r2r_2d(d1,d2,f->data(),h->data(), FFTW_REDFT00, FFTW_REDFT00, FFTW_ESTIMATE);
}

void Grid::populate_iid() {
    srand(time(NULL));
    for(int i=0; i < d1*d2; i++ )
        h->at(i) = ((double)rand()) / RAND_MAX;
}

void Grid::populate_coscos(double fact1, double fact2) {
    for(int i=0; i < d1*d2; i++ )
        h->at(i) = cos( fact1 * ( i % d1 ) ) * cos( fact2 * ( i / d1 ) );
}

void Grid::populate_square(int side) {
    for( int i= ( d1 - side ) / 2; i < ( d1 + side) / 2 + 1; i++ ) 
        for( int j= ( d2 - side ) / 2; j < ( d2 + side) / 2 + 1; j++ )
            h->at( i * d1 + j ) = 0.5;
}

void Grid::fourier_transform(double reverse) {
    if(reverse) fftw_execute( reverse_plan );
    else fftw_execute( direct_plan );
}

void Grid::multiply_fft(double gamma) {
    f->at(0) = 0;
    for(int i=1; i<d1*d2; i++)
        f->at(i) = f->at(i) * pow( d(i,0), (gamma-1)/2 );
}

void Grid::print_data(const char* filename, double fourier) const {
    ofstream out(filename);
    for(int i=0; i<d1; i++)
        for(int j=0; j<d2; j++)
            out<<i<<'\t'<<j<<'\t'<<operator()(i,j,fourier)<<'\n';
    out.close();
}

double Grid::operator()(const int i) const {
    return h->at( i );
}

double Grid::operator()(const int i, const int j, double fourier) const {
    if( fourier ) return f->at( d1 * i + j );
    return h->at( d1 * i + j );
}

const double Grid::d(int i1, int i2) const{
    return d( i1 / d1, i1 % d1, i2 / d1, i2 % d2 );
}

const double Grid::d(int x1, int y1, int x2, int y2) const{
    return sqrt( pow( x1 - x2, 2 ) + pow( y1 - y2, 2 ) );
}

vector< tuple<double,double, double> >* Grid::corr_func(double resol) const {
    vector< double > sum( floor( sqrt( d1 * d1 + d2 * d2 ) / resol ) + 1, 0 );
    vector< double > sum2( sum.size(), 0 );
    vector< int    > count( sum.size(), 0 );
    int d_;
    for( int i = 0; i < d1*d2; i++ ) {
        for( int j = 0; j < d1*d2; j++ ) {
            d_ = d( i, j );
            sum [ d_ ] += operator()(i) * operator()(j);
            sum2[ d_ ] += operator()(i) * operator()(j) * operator()(i) * operator()(j);
            count[ d_ ] ++;
        }
        if( i % d1 == 0 ) cout<<i/d1<<endl;
    }
    auto r = new vector< tuple< double, double, double > >( sum.size() );
    for( int i = 0; i < sum.size(); i++ ) 
        r->at(i) = make_tuple( resol * i, sum[i] / count[i], sqrt( ( sum2[i] - pow( sum[i], 2 ) / count[i] ) / ( count[i] - 1 ) ) );
    return r;
}

void Grid::print_corr(const char* filename) const {
    ofstream out(filename);
    auto cf = corr_func(1);
    for( auto d : *cf) {
        out<<get<0>(d)<<'\t'<<get<1>(d)<<'\t'<<get<2>(d)<<'\n';
    }
    out.close();
}