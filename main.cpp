#include <cstdlib>
#include <iostream>
#include "Grid.h"
#include "CorrFuncCalcolator.h"
#include "FourierCoupledGrids.h"
#include "GridFiller.h"
#include "date.h"
#include <chrono>
#include <fstream>
#include <filesystem>

using namespace std;

void populate_defects(Grid<int>& g, FourierCoupledGrids& fcg, double gamma, int howmuch) {
    GridFiller::clean(g);

    GridFiller::iid( fcg.h );
    fcg.fourier_transform();
    fcg.multiply_fft_new( gamma );
    fcg.fourier_transform(true);

    GridFiller::ranked_insertion( g, fcg.h, howmuch);
}

void run_replies(int side, double defects_frac, double gamma, int n_replies, int corr_range, string path) {
    corr_range = min( side/2, corr_range );
    vector< pair<int, double> > corr_func_avg_h ( corr_range, make_pair(0,0) );
    vector< pair<int, double> > corr_func_avg_d ( corr_range, make_pair(0,0) );

    FourierCoupledGrids fcg(side);
    Grid<int> g(side);
    CorrFunc::Expospaced<double> CF_H( &fcg.h, 1.4, corr_range, 1 );
    CorrFunc::Expospaced<int>    CF_D( &g    , 1.4, corr_range, 1 );

    for( int j=0; j < corr_range; j++ ) {
        corr_func_avg_h[j].first = j;
        corr_func_avg_d[j].first = j;
    }

    GridFiller::iid( fcg.h );

    for( int i=0; i < n_replies; i++ ) {

        populate_defects( g, fcg, gamma, defects_frac * side * side );

        if( corr_range > 0 ) {
            vector< CorrFunc::Datapoint >* cfh = CF_H.compute_corr_function ( );
            vector< CorrFunc::Datapoint >* cfd = CF_D.compute_corr_function ( );

            for( int j=0; j < min( corr_range, min( (int)cfh->size(), (int)cfd->size() ) ); j++ ) {
                corr_func_avg_h[j].second = corr_func_avg_h[j].second + cfh->at(j).value;
                corr_func_avg_d[j].second = corr_func_avg_d[j].second + cfd->at(j).value;
            }
            
            cfh->clear();
            delete cfh;
            cfd->clear();
            delete cfd;
        }
        
        if( i % ( n_replies/100 + 1 ) == 0 ) cout<<i<<endl;
    }

    if( !filesystem::exists( path ) )
        filesystem::create_directory( path );

    ofstream out_details( path + "/details.txt");
    out_details<<"{\n\"side\":\t"<<side<<",\n\"defects_frac\":\t"<<defects_frac<<
               ",\n\"gamma\":\t"<<gamma<<",\n\"replies\":"<<n_replies<<
               ",\n\"corr_range\":\t"<<corr_range<<"\n}"<<endl;

    if( corr_range > 0 ) {
        ofstream out_corr( path + "/corr_func_avg_h.txt");
        for( int j=0; j < corr_range; j++ ) {
            corr_func_avg_h[j].second /= n_replies;
            out_corr<< corr_func_avg_h[j].first << '\t' << corr_func_avg_h[j].second << '\n';
        }

        ofstream out2_corr( path + "/corr_func_avg_d.txt");
        for( int j=0; j < corr_range; j++ ) {
            corr_func_avg_d[j].second /= n_replies;
            out2_corr<< corr_func_avg_d[j].first << '\t' << corr_func_avg_d[j].second << '\n';
        }
    }

}

int main() {

    int side = 1000;
    int corr_range = side/2; // 30
    double defects_frac = 0.3;
    double gamma = 0.4;
    int n_replies = 1;
    string path;
    int progress = 0;
    bool ask = false;

    do {
        path = date::format("%Y%m%d", chrono::system_clock::now()) + "_" + to_string( progress++ );
    } while( filesystem::exists(path) );

    if( ask ) {
        cout<<"Choose side (suggested "<<side<<"): ";
        cin>>side;
        cout<<"Choose defects frac (suggested "<<defects_frac<<"): ";
        cin>>defects_frac;
        cout<<"Choose gamma (suggested "<<gamma<<"): ";
        cin>>gamma;
        cout<<"Choose number of replies (suggested "<<n_replies<<"): ";
        cin>>n_replies;
        cout<<"Choose correlation range (suggested "<<corr_range<<"): ";
        cin>>corr_range;
    }

    cout<<"\aStarting..."<<endl;

    run_replies( side, defects_frac, gamma, n_replies, corr_range, path );

    cout<<"\a\a\aDone!"<<endl;
    system("PAUSE");
}