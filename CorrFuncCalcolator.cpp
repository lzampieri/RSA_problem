#include "CorrFuncCalcolator.h"

using namespace std;
using namespace CorrFunc;

template<class T>
void BaseClass<T>::_thread_postman( ) {
    works_i_mutex.lock();
    while( works_i < works->size() ) {
        Work* todo = works->at( works_i );
        works_i += 1;
        works_i_mutex.unlock();

        _thread_worker( todo );

        works_i_mutex.lock();
    }
    works_i_mutex.unlock();
}

template<class T>
void BaseClass<T>::_thread_worker( Work* work ) {
    double sum = 0, sum2 = 0;
    int count = 0;
    int v1 = work->v1, v2 = work->v2;

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

    _thread_update_rawdata(
        RawDatapoint(sum, sum2, count),
        work->raw_data_i
    );
}

template<class T>
void BaseClass<T>::_thread_update_rawdata( RawDatapoint rdp, int raw_data_i ) {
    raw_data_mutex.lock();
    raw_data->at(raw_data_i) += rdp;
    raw_data_mutex.unlock();
}

template<class T>
BaseClass<T>::BaseClass( const Grid<T>* grid ) : grid(grid) {
    works = new vector< Work* >();
    is = new vector< int >();
    raw_data = new vector< RawDatapoint >();
    final_data = new vector< Datapoint >();
}

template<class T>
void BaseClass<T>::auto_populate_works(map<int,int>& to_add, int max_v) {
    if( max_v == -1 ) max_v = grid->d1;
    for(int v1 = 0; v1 < max_v; v1++ ) {
        for(int v2 = v1; v2 < max_v; v2++) {
            auto i = to_add.find( grid->d(0,0,v1,v2) );
            if( i != to_add.end() ) {
                this->works->push_back( new Work( v1, v2, i->second ) );
            }
        }
    }
    cout<<"Number of works to be carried out: "<<this->works->size()<<endl;
}

template<class T>
vector< Datapoint >* BaseClass<T>::compute_corr_function() {
    vector<thread*> threads;
    threads.reserve( CORRFUNC_MAX_THREADS );

    // Clean raw_data
    for( int i = 0; i < raw_data->size(); i++ ) {
        raw_data->at(i).sum = 0;
        raw_data->at(i).sum2 = 0;
        raw_data->at(i).count = 0;
    }

    works_i = 0;

    // Launch threads
    for( int i = 0; i < CORRFUNC_MAX_THREADS; i++ ) {
        threads.push_back( 
                new thread( _thread_postman, this )
            );
    }

    // Joint threads
    for( int i = 0; i < threads.size(); i++ ) {
        threads[i]->join();
    }
    threads.clear();

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
void BaseClass<T>::print_corr( const char* filename ) {
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

template<class T>
Equispaced<T>::Equispaced(const Grid<T>* grid, int start, int end, int step) : BaseClass<T>(grid) {
    map<int,int> themap;

    for( int i = start; i < end; i+=step ) {
        this->is->push_back( i );
        this->final_data->push_back( Datapoint( 0, 0 ) );
        this->raw_data->push_back( RawDatapoint( 0, 0, 0 ) );
        themap[i] = this->final_data->size() - 1 ;
    }

    this->auto_populate_works(themap, end);
}

template<class T>
Expospaced<T>::Expospaced(const Grid<T>* grid, double base, int maxval, int exp_step) : BaseClass<T>(grid) {
    map<int,int> themap;

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
        this->is->push_back( val );
        this->final_data->push_back( Datapoint( 0, 0 ) );
        this->raw_data->push_back( RawDatapoint( 0, 0, 0 ) );
        themap[val] = this->final_data->size() - 1 ;
    } while( update_val() );

    this->auto_populate_works(themap, maxval);
}