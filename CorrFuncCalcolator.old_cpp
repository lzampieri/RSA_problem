#include "CorrFuncCalcolator.h"

using namespace std;
using namespace CorrFunc;

template<class T>
void Calculator<T>::run_works( ) {
    // Questa struttura un po' strana è eredità del fatto che una volta questa cosa era a thread
    while( works_i < works->size() ) {
        run_work( works->at( works_i ) );
        works_i += 1;
    }
}

template<class T>
void Calculator<T>::run_work( Work* work ) {
    double sum = 0, sum2 = 0;
    int count = 0;
    GridSite v = work->v;
    GridSite u = GridSite( v.Y, v.X, v );
    GridSite vR = GridSite( v.X, -v.Y, v );
    GridSite uR = GridSite( v.Y, -v.X, v );

    GridSite xy1( v ), xy2( v ), xy3( v ), xy4( v ), xy5( v ); // Random initialization to allocate space
    double val1, val2, val3, val4, val5;
    for( int i = 0; i < grid->d1 * grid->d2; i++ ) {
        GridSite xy1 = grid->_xy( i );
        GridSite xy2 = xy1 + v;
        GridSite xy3 = xy1 + u;
        GridSite xy4 = xy1 + vR;
        GridSite xy5 = xy1 - uR;

        val1 = grid->operator[]( xy1 );
        val2 = grid->operator[]( xy2 );
        val3 = grid->operator[]( xy3 );
        val4 = grid->operator[]( xy4 );
        val5 = grid->operator[]( xy5 );
        
        sum += val1*val2;
        sum2+= val1*val2*val1*val2;
        count++;
        
        if( v != u ) {
            sum += val1*val3;
            sum2+= val1*val3*val1*val3;
            count++;
        }
        
        sum += val1*val4;
        sum2+= val1*val4*val1*val4;
        count++;
        
        if( v != u ) {
            sum += val1*val5;
            sum2+= val1*val5*val1*val5;
            count++;
        }
    }

    update_data(
        RawDatapoint(sum, sum2, count),
        work->raw_data_i
    );
}

template<class T>
void Calculator<T>::update_data( RawDatapoint rdp, int raw_data_i ) {
    raw_data->at(raw_data_i) += rdp;
}

template<class T>
Calculator<T>::Calculator( Model* model, const Grid<T>* grid ) : grid(grid) {
    works = new vector< Work* >();
    is = new vector< int >();
    map<int,int> themap;
    int max_dist = grid->d1;

    for( int i=0; i < model->is.size(); i++ ) {
        is->push_back( model->is[i] );
        themap[ model->is[i] ] = i;
        max_dist = max( max_dist, model->is[i] );
    }

    for(int v1 = 0; v1 < max_dist; v1++ ) {
        for(int v2 = v1; v2 < max_dist; v2++) {
            auto i = themap.find( grid->d(0,0,v1,v2) );
            if( i != themap.end() ) {
                works->push_back( new Work( GridSite( v1, v2, *grid ), i->second ) );
            }
        }
    }

    raw_data = new vector< RawDatapoint >( model->is.size(), RawDatapoint( 0, 0, 0 ) );
    final_data = new vector< Datapoint >( model->is.size(), Datapoint( 0, 0 ) );
}

template<class T>
Calculator<T>::~Calculator() {
    delete works;
    delete is;
    delete raw_data;
    delete final_data;
}

template<class T>
const vector< Datapoint >* Calculator<T>::compute_corr_function() {

    // Clean raw_data
    for( int i = 0; i < raw_data->size(); i++ ) {
        raw_data->at(i).sum = 0;
        raw_data->at(i).sum2 = 0;
        raw_data->at(i).count = 0;
    }

    works_i = 0;

    // Launch analysis
    run_works();

    // Compute real data
    for( int i = 0; i < raw_data->size(); i++ ) {
        final_data->at(i).value = raw_data->at(i).sum / raw_data->at(i).count;
        final_data->at(i).std   = sqrt(
                            ( raw_data->at(i).sum2 - pow( raw_data->at(i).sum, 2 ) / raw_data->at(i).count )
                            /
                            ( raw_data->at(i).count - 1 )
                        );
    }

    return final_data;
}

template<class T>
void Calculator<T>::print_corr( const char* filename ) {
    ostream* out;
    if( strlen( filename ) > 0 )
        out = new ofstream(filename);
    else
        out = &cout;

    compute_corr_function();

    for( int i = 0; i < final_data->size(); i++ ) {
        (*out)<< is->at(i) <<'\t'<< final_data->at(i).value <<'\t'<< final_data->at(i).std <<'\n';
    }

    delete out;
}

Equispaced::Equispaced(int start, int end, int step) {
    for( int i = start; i < end; i+=step ) {
        is.push_back( i );
    }
    keyname = "Equispaced " + to_string( end );
}

Expospaced::Expospaced(double base, int maxval, int exp_step) {
    int exp = 0;
    int val = 1;
    auto update_val = [&exp, &val, maxval, exp_step, base] () {
        int old_val = val;
        do {
            exp+=exp_step;
            val = pow(base,exp);
        } while( old_val == val );
        return val < maxval;
    };

    do {
        is.push_back( val );
    } while( update_val() );
    keyname = "Expospaced " + to_string( maxval );
}