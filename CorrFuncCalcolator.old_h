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
#include <string>
#include "Grid.h"

#include <iostream>

namespace CorrFunc {

    struct Work {
        GridSite v;
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

    class Model {
    public:
        std::vector< int > is;
        std::string keyname;
    };

    template<class T>
    class Calculator {

    private:
        int works_i;

        void run_works( );
        void run_work( Work* work );
        void update_data( RawDatapoint rdp, int raw_data_i );

        const Grid<T>* grid;

    protected:
        std::vector< Work* >* works;
        std::vector< RawDatapoint >* raw_data;
        std::vector< Datapoint >* final_data;

    public:
        Calculator( Model* model, const Grid<T>* grid );
        ~Calculator();
        
        std::vector< int >* is;
        const std::vector< Datapoint >* compute_corr_function();
        void print_corr(const char* filename);

    };

    class Equispaced : public Model {
    public:
        Equispaced(int start, int end, int step);
    };

    class Expospaced : public Model {
    public:
        Expospaced(double base, int maxval, int exp_step);
    };

    template class Calculator<double>;
    template class Calculator<int>;  

}

#endif