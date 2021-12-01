#include "CorrFuncCalcolator.h"
#include "Grid.h"

void CorrFuncCalcolator::_thread_worker(
        const Grid* grid, const int v1, const int v2,
        vector< CorrFunc_datapoint >* results, mutex* m ) {

    double sum = 0, sum2 = 0;
    int count = 0;

    pair<int, int> xy1, xy2, xy3;
    double val1, val2, val3;
    for( int i = 0; i < grid->d1 * grid->d2; i++ ) {
        xy1 = grid->_xy( i );
        xy2 = make_pair( xy1.first + v1, xy1.second + v2 );
        xy3 = make_pair( xy1.first + v2, xy1.second + v1 );

        val1 = grid->operator()( xy1 );
        val2 = grid->operator()( xy2 );
        val3 = grid->operator()( xy3 );
        
        sum += val1*val2;
        sum2+= val1*val2*val1*val2;
        count++;
        
        sum += val1*val3;
        sum2+= val1*val3*val1*val3;
        count++;
    }

    m->lock();
    results->push_back( make_tuple(
        grid->d( 0, 0, v1, v2 ),
        sum / count,
        sqrt( ( sum2 - pow( sum, 2 ) / count ) / ( count - 1 ) )
    ) );
    m->unlock();
}

vector< tuple<double,double, double> >* CorrFuncCalcolator::compute_corr_function(const Grid* grid) {
    queue<thread*> threads;
    mutex* m = new mutex();
    vector< tuple<double, double, double> >* results = new vector< tuple<double, double, double> >();

    for( int v1 = 0; v1 < grid->d1 / 2; v1++ ) {
        for( int v2 = v1; v2 < grid->d2 / 2; v2++ ) {

            threads.push( 
                new thread( _thread_worker, grid, v1, v2, results, m )
            );

        }
    }

    for( int i = 0; i < threads.size(); i++ ) {
        threads.front()->join();
        threads.pop();
    }

    return results;
}

void CorrFuncCalcolator::print_corr(const Grid* grid, const char* filename) {
    ofstream out(filename);
    auto cf = compute_corr_function(grid);
    for( auto d : *cf) {
        out<<get<0>(d)<<'\t'<<get<1>(d)<<'\t'<<get<2>(d)<<'\n';
    }
    out.close();
}
void CorrFuncCalcolator::display_corr(const Grid* grid) {
    auto cf = compute_corr_function(grid);
    for( auto d : *cf) {
        cout<<get<0>(d)<<'\t'<<get<1>(d)<<'\t'<<get<2>(d)<<'\n';
    }
}