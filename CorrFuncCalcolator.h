#ifndef CORRFUNCCALCOLATOR_H
#define CORRFUNCCALCOLATOR_H

#include <cstdlib>
#include <vector>
#include <map>
#include <cmath>
#include <fstream>
#include <algorithm>
#include <climits>
#include <cstring>
#include "Grid.h"

#include <thread>
#include <mutex>

#define CORRFUNC_MAX_THREADS 30

#include <iostream>

namespace CorrFunc {

    struct Work {
        int v1;
        int v2;
        int raw_data_i;
    };

    struct RawDatapoint {
        double sum;
        double sum2;
        double count;

        RawDatapoint& operator+=(const RawDatapoint& rhs) {
            sum += rhs.sum;
            sum2 += rhs.sum2;
            count += rhs.count;
            return *this;
        }
    };

    struct Datapoint {
        double value;
        double std;
    };

    template<class T>
    class BaseClass {

    private:
        std::mutex works_i_mutex;
        int works_i;
        std::mutex raw_data_mutex;

        void _thread_postman( );
        void _thread_worker( Work* work );
        void _thread_update_rawdata( RawDatapoint rdp, int raw_data_i );

        const Grid<T>* grid;

    protected:
        std::vector< Work* >* works;
        std::vector< RawDatapoint >* raw_data;
        std::vector< Datapoint >* final_data;

        BaseClass( const Grid<T>* grid );

        void auto_populate_works(std::map<int,int>& to_add,int max_v=-1);

    public:
        ~BaseClass();
        
        std::vector< int >* is;
        const std::vector< Datapoint >* compute_corr_function();
        void print_corr(const char* filename);

    };

    template<class T>
    class Equispaced : public BaseClass<T> {
    public:
        Equispaced(const Grid<T>* grid, int start, int end, int step);
    };

    template<class T>
    class Expospaced : public BaseClass<T> {
    public:
        Expospaced(const Grid<T>* grid, double base, int maxval, int exp_step);
    };

    template class BaseClass<double>;
    template class BaseClass<int>;    
    template class Equispaced<double>;
    template class Equispaced<int>;    
    template class Expospaced<double>;
    template class Expospaced<int>;

}

#endif