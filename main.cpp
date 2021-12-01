#include <cstdlib>
#include <iostream>
#include "Grid.h"
#include "CorrFuncCalcolator.h"
#include <fstream>

using namespace std;

int main() {
    Grid g(100);
    g.populate_iid();

    CorrFuncCalcolator::print_corr( &g, "../data/cf1.txt" );

    // // g.print_data("../data/real.txt");

    g.fourier_transform();

    // // g.print_data("../data/fft.txt", true);

    g.multiply_fft( 0.02 );

    // // g.print_data("../data/fft2.txt", true);

    g.fourier_transform(true);

    // // g.print_data("../data/real2.txt");
    
    CorrFuncCalcolator::print_corr( &g, "../data/cf2.txt" );

    system("PAUSE");
}