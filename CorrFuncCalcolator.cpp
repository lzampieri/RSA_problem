#include "CorrFuncCalcolator.h"
#include "Grid.h"

using namespace std;

void CorrFuncCalcolator::_thread_postman(
        const Grid* grid, queue< CorrFunc_Work >* todolist,
        vector< CorrFunc_RawDatapoint >* results, mutex* m ) {
    m->lock();
    while( !todolist->empty() ) {
        CorrFunc_Work todo = todolist->front();
        todolist->pop();
        m->unlock();

        CorrFunc_RawDatapoint r = _thread_worker( grid, todo.first, todo.second );

        m->lock();
        results->push_back( r );
    }

    m->unlock();
}

CorrFunc_RawDatapoint CorrFuncCalcolator::_thread_worker(
        const Grid* grid, const int v1, const int v2 ) {

    double sum = 0, sum2 = 0;
    int count = 0;

    pair<int, int> xy1, xy2, xy3, xy4, xy5;
    double val1, val2, val3, val4, val5;
    for( int i = 0; i < grid->d1 * grid->d2; i++ ) {
        xy1 = grid->_xy( i );
        xy2 = make_pair( xy1.first + v1, xy1.second + v2 );
        xy3 = make_pair( xy1.first + v2, xy1.second + v1 );
        xy4 = make_pair( xy1.first + v1, xy1.second - v2 );
        xy5 = make_pair( xy1.first + v2, xy1.second - v1 );

        val1 = grid->operator()( xy1 );
        val2 = grid->operator()( xy2 );
        val3 = grid->operator()( xy3 );
        val4 = grid->operator()( xy4 );
        val5 = grid->operator()( xy5 );
        
        sum += val1*val2;
        sum2+= val1*val2*val1*val2;
        count++;
        
        if( v1 != v2 ) {
            sum += val1*val3;
            sum2+= val1*val3*val1*val3;
            count++;
        }
        
        sum += val1*val4;
        sum2+= val1*val4*val1*val4;
        count++;
        
        if( v1 != v2 ) {
            sum += val1*val5;
            sum2+= val1*val5*val1*val5;
            count++;
        }
    }

    return make_tuple(
        grid->d( 0, 0, v1, v2 ),
        count,
        sum,
        sum2
    );
}

vector< CorrFunc_Datapoint >* CorrFuncCalcolator::compute_corr_function(const Grid* grid, int max_range) {
    queue<thread*> threads;
    queue<CorrFunc_Work> works;
    mutex* m = new mutex();
    vector< CorrFunc_RawDatapoint >* raw_results = new vector< CorrFunc_RawDatapoint >();

    for( int v1 = 0; v1 < max_range; v1++ ) {
        for( int v2 = v1; v2 < max_range; v2++ ) {
            works.push( make_pair( v1, v2) );
        }
    }

    for( int i = 0; i < CORRFUNC_MAX_THREADS; i++ ) {
        threads.push( 
                new thread( _thread_postman, grid, &works, raw_results, m )
            );
    }

    for( int i = 0; i < threads.size(); i++ ) {
        threads.front()->join();
        threads.pop();
    }

    sort(
        raw_results->begin(),
        raw_results->end(),
        [](const CorrFunc_RawDatapoint p1, const CorrFunc_RawDatapoint p2) { return get<0>(p1) < get<0>(p2); } );

    vector< CorrFunc_Datapoint >* results = new vector< CorrFunc_Datapoint >();
    int d = 0;
    double sum = 0, sum2 = 0;
    int count = 0;
    for( CorrFunc_RawDatapoint dp : *raw_results ) {
        if( (int)get<0>(dp) > d ) {
            if( count > 0 ) {
                results->push_back( make_tuple (
                    d,
                    sum / count,
                    sqrt( ( sum2 - pow( sum, 2 ) / count ) / ( count - 1 ) )
                ));
                sum = 0;
                sum2 = 0;
                count = 0;
                d = (int)get<0>(dp);
            }
        }

        count+= get<1>(dp);
        sum += get<2>(dp);
        sum2+= get<3>(dp);
    }

    raw_results->clear();
    delete raw_results;
    delete m;

    return results;
}

void CorrFuncCalcolator::print_corr(const Grid* grid, const char* filename, int max_range) {
    ofstream out(filename);
    vector< CorrFunc_Datapoint >* cf = compute_corr_function(grid, max_range);
    for( CorrFunc_Datapoint d : *cf) {
        out<<get<0>(d)<<'\t'<<get<1>(d)<<'\t'<<get<2>(d)<<'\n';
    }

    cf->clear();
    delete cf;
    out.close();
}
void CorrFuncCalcolator::display_corr(const Grid* grid, int max_range) {
    auto cf = compute_corr_function(grid, max_range);
    for( auto d : *cf) {
        cout<<get<0>(d)<<'\t'<<get<1>(d)<<'\t'<<get<2>(d)<<'\n';
    }

    cf->clear();
    delete cf;
}