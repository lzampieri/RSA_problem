#include <cstdlib>
#include <iostream>
#include "Grid.h"
#include "CorrFuncCalcolator.h"
#include "FourierCoupledGrids.h"
#include "GridFiller.h"
#include <fstream>

using namespace std;

void populate_defects(Grid<int>& g, FourierCoupledGrids& fcg, double gamma, int howmuch) {
    GridFiller::clean(g);

    GridFiller::iid( fcg.h );
    fcg.fourier_transform();
    fcg.multiply_fft_new( gamma );
    fcg.fourier_transform(true);

    GridFiller::ranked_insertion( g, fcg.h, howmuch);
}

// #include "fftw3.h"

int main() {

    int side = 1024;
    int corr_range = 30;
    double defects_frac = 0.3;
    double gamma = 0.4;
    int N = 10;

    // ofstream out1("../data/powerspectra_test.txt");

    // vector<double> f1(side,0), f2(side,0), h1(side,0), h2(side,0);
    // vector<double> h1_(side,0), h2_(side,0), f1_(side,0), f2_(side,0);
    // fftw_plan plan1 = fftw_plan_r2r_1d( side, f1.data(), h1.data(), FFTW_REDFT00, FFTW_ESTIMATE),
    //           plan2 = fftw_plan_r2r_1d( side, f2.data(), h2.data(), FFTW_REDFT00, FFTW_ESTIMATE),
    //           plan1_rev = fftw_plan_r2r_1d( side, h1_.data(), f1_.data(), FFTW_REDFT00, FFTW_ESTIMATE),
    //           plan2_rev = fftw_plan_r2r_1d( side, h2_.data(), f2_.data(), FFTW_REDFT00, FFTW_ESTIMATE);

    // int d = 1;
    // bool radq = false;

    // for(double q = 0; q < side; q++ ) {

    //     // Populate f with the expected power law
    //     double Sold = pow( q, (gamma-1) );
    //     double q1 = radq ? sqrt( Sold ) : Sold;
    //     double beta = ( gamma - 1 ) / 2;
    //     double pref = tgamma( beta + 1 );
    //     double qrn = M_PI * q / ( side - 1 );
    //     double S = pref * pow( qrn / 2 , beta ) * cyl_bessel_k( abs(beta), qrn );
    //     double q2 = radq ? sqrt( S ) : S;
    //     f1[q] = q1;
    //     f2[q] = q2;

    //     // Populate h_ with the expected corr. function
    //     h1_[q] = pow( q, -gamma );
    //     h2_[q] = pow( 1 + q*q, -gamma / 2 );
    // }
    // h1_[0] = 1;
    // f1[0] = 0;
    // f2[0] = 0;
    // fftw_execute( plan1 );
    // fftw_execute( plan2 );
    // fftw_execute( plan1_rev );
    // fftw_execute( plan2_rev );
    // for(double q = 0; q < side; q++ ) {
    //     out1<<q<<'\t'<<f1[q]<<'\t'<<f2[q]<<'\t'<<h1[q]<<'\t'<<h2[q]<<'\t'<<f1_[q]<<'\t'<<f2_[q]<<'\t'<<h1_[q]<<'\t'<<h2_[q]<<'\n';
    // }
    // return 0;

    // FourierCoupledGrids fcg(side);
    // GridFiller::iid( fcg.h );
    // fcg.h.print_data("../data/h1.txt");
    // CorrFuncCalcolator<double>::print_corr( &fcg.h, "../data/cf1.txt", corr_range );
    // fcg.fourier_transform();
    // fcg.multiply_fft_new( gamma );
    // fcg.fourier_transform(true);
    // fcg.h.print_data("../data/h2.txt");
    // CorrFuncCalcolator<double>::print_corr( &fcg.h, "../data/cf2.txt", corr_range );

    // return 0;

    corr_range = min( side, corr_range );
    vector< pair<int, double> > corr_func_avg_h ( corr_range, make_pair(0,0) );
    vector< pair<int, double> > corr_func_avg_d ( corr_range, make_pair(0,0) );

    FourierCoupledGrids fcg(side);
    Grid<int> g(side);

    for( int j=0; j < corr_range; j++ ) {
        corr_func_avg_h[j].first = j;
        corr_func_avg_d[j].first = j;
    }

    GridFiller::iid( fcg.h );

    for( int i=0; i < N; i++ ) {

        populate_defects( g, fcg, gamma, defects_frac * side * side );
        vector< CorrFunc_Datapoint >* cfh = CorrFuncCalcolator<double>::compute_corr_function ( &(fcg.h), corr_range );
        vector< CorrFunc_Datapoint >* cfd = CorrFuncCalcolator<int>::compute_corr_function ( &g, corr_range );

        // fcg.h.print_data("../data/fcg_h.txt");
        // g.print_data("../data/g.txt");

        for( int j=0; j < min( corr_range, min( (int)cfh->size(), (int)cfd->size() ) ); j++ ) {
            corr_func_avg_h[j].second = corr_func_avg_h[j].second + get<1>(cfh->at(j));
            corr_func_avg_d[j].second = corr_func_avg_d[j].second + get<1>(cfd->at(j));
        }

        cfh->clear();
        delete cfh;
        cfd->clear();
        delete cfd;
        
        if( i % ( N/100 + 1 ) == 0 ) cout<<i<<endl;
    }

    ofstream out("../data/corr_func_avg_h.txt");
    for( int j=0; j < corr_range; j++ ) {
        corr_func_avg_h[j].second /= N;
        out<< corr_func_avg_h[j].first << '\t' << corr_func_avg_h[j].second << '\n';
    }

    ofstream out2("../data/corr_func_avg_d.txt");
    for( int j=0; j < corr_range; j++ ) {
        corr_func_avg_d[j].second /= N;
        out2<< corr_func_avg_d[j].first << '\t' << corr_func_avg_d[j].second << '\n';
    }

    // bool corr_func = true;
    // FourierCoupledGrids g(side);
    // cout<<"Created"<<endl;

    // GridFiller::iid( g.h );
    // cout<<"Populated"<<endl;

    // if( corr_func ) {
    //     CorrFuncCalcolator<double>::print_corr( &g.h, "../data/cf1.txt", min(side,corr_range) );
    //     cout<<"Corr Func calcolated"<<endl;
    // }

    // g.h.print_data("../data/real.txt");

    // g.fourier_transform();

    // g.multiply_fft_new( 0.4 );

    // g.fourier_transform(true);
    // cout<<"Correlated"<<endl;

    // g.h.print_data("../data/real2.txt");

    // if( corr_func ) {    
    //     CorrFuncCalcolator<double>::print_corr( &g.h, "../data/cf2.txt", min(side,corr_range) );
    //     cout<<"Corr Func calcolated"<<endl;
    // }

    // Grid<int> g2( side );
    // cout<<"Defects grid created"<<endl;

    // GridFiller::montecarlo( g2, g.h, defects_frac * side * side );
    // cout<<"Grid filled"<<endl;

    // g2.print_data("../data/defects.txt");

    // if( corr_func ) {    
    //     CorrFuncCalcolator<int>::print_corr( &g2, "../data/cf_defects.txt", min(side,corr_range) );
    //     cout<<"Corr Func calcolated"<<endl;
    // }

    system("PAUSE");
}