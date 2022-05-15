#include "NewCFH.h"

using namespace std;
using namespace NewCFH;

// Constructor
Calculator::Calculator( NewCF::Model* model, Grid<double>* grid ) :
    model( model ), grid( grid ) {
    values = new vector<double>( model->items.size() );
}

// Worker
vector< double >* Calculator::calculate() {
    for( int i = 0; i < values->size(); i++ ) {
        values->at( i ) = 0;
    }

    for( int i_i = 0; i_i < values->size(); i_i++ ) {
        for( int i_v = 0; i_v < model->items[i_i].second.size(); i_v++ ) {
            int count = 0;
            for( int i_g = 0; i_g < grid->imax(); i_g ++ ) {
                GridSite here = GridSite( i_g, *grid );
                int new_x = here.X + model->items.at(i_i).second[i_v].first;
                int new_y = here.Y + model->items.at(i_i).second[i_v].second;
                if( ( new_x > here.d1 ) || ( new_x < 0 ) || ( new_y > here.d2 ) || ( new_y < 0 ) )
                    continue;
                count++;
                GridSite there = GridSite( new_x, new_y, *grid );
                // GridSite there = GridSite( i_g, *grid )
                //     + GridSite(
                //         model->items.at(i_i).second[i_v].first,
                //         model->items.at(i_i).second[i_v].second,
                //         *grid
                //         );
                values->at( i_i ) += grid->operator[]( here ) * grid->operator[](  there );
            }
            values->at( i_i ) /= count;
            // values->at( i_i ) /= grid->imax();
            // cout<<model->items[i_i].first<<'\t'<<model->items.at(i_i).second[i_v].first<<'\t'<<model->items.at(i_i).second[i_v].second<<endl;
        }
    }

    for( int i = 0; i < values->size(); i++ ) {
        values->at( i ) /= model->items[i].second.size();
    }

    return values;
}