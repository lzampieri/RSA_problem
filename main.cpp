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

    GridFiller::montecarlo( g, fcg.h, howmuch);
}

int main() {
    int side = 200;
    int corr_range = 30;
    double defects_frac = 0.2;
    double gamma = 0.4;
    int N = 2000;

    ofstream out("../data/cf_defects.txt");

    corr_range = min( side, corr_range );
    vector< pair<int, double> > corr_func_avg ( corr_range, make_pair(0,0) );

    FourierCoupledGrids fcg(side);
    Grid<int> g(side);

    for( int j=0; j < corr_range; j++ ) corr_func_avg[j].first = j;

    GridFiller::iid( fcg.h );
    CorrFuncCalcolator<double>::print_corr( &fcg.h, "../data/cf1.txt", corr_range );

    for( int i=0; i < N; i++ ) {

        populate_defects( g, fcg, gamma, defects_frac * side * side );
        vector< CorrFunc_Datapoint >* cfd = CorrFuncCalcolator<int>::compute_corr_function ( &g, corr_range );

        if( i == 0 )
            CorrFuncCalcolator<double>::print_corr( &fcg.h, "../data/cf2.txt", corr_range );

        for( int j=0; j < min( corr_range, (int)cfd->size() ); j++ ) {
            corr_func_avg[j].second = corr_func_avg[j].second + get<1>(cfd->at(j));
        }

        cfd->clear();
        delete cfd;
        
        if( i % 10 == 0 ) cout<<i<<endl;
    }

    for( int j=0; j < corr_range; j++ ) {
        corr_func_avg[j].second /= N;
        out<< corr_func_avg[j].first << '\t' << corr_func_avg[j].second << '\n';
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