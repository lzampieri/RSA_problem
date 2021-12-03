#include <cstdlib>
#include <iostream>
#include "Grid.h"
#include "CorrFuncCalcolator.h"
#include <fstream>

using namespace std;

int main() {
    Grid g(500);
    cout<<"Created"<<endl;

    GridFiller::iid(g);
    cout<<"Populated"<<endl;

    CorrFuncCalcolator::print_corr( &g, "../data/cf1.txt", 50 );
    cout<<"Corr Func calcolated"<<endl;

    // g.print_data("../data/real.txt");

    g.fourier_transform();

    // g.print_data("../data/fft.txt", true);

    g.multiply_fft_new( 0.4 );

    // g.print_data("../data/fft2.txt", true);

    g.fourier_transform(true);
    cout<<"Correlated"<<endl;

    // g.print_data("../data/real2.txt");
    
    CorrFuncCalcolator::print_corr( &g, "../data/cf2.txt", 50 );
    cout<<"Corr Func calcolated"<<endl;

    system("PAUSE");
}