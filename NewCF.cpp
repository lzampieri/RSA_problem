#include "NewCF.h"

using namespace std;
using namespace NewCF;

// Constructor
Calculator::Calculator( Model* model, Grid<int>* grid ) :
    model( model ), grid( grid ) {
    values = new vector<double>( model->items.size() );
}

// Worker
vector< double >* Calculator::calculate() {
    for( int i = 0; i < values->size(); i++ ) {
        values->at( i ) = 0;
    }
    double gshere, gsthere;
    for( int i_i = 0; i_i < values->size(); i_i++ ) {
        for( int i_v = 0; i_v < model->items[i_i].second.size(); i_v++ ) {
            int count = 0;
            for( int i_g = 0; i_g < grid->imax(); i_g ++ ) {
                GridSite here = GridSite( i_g, *grid );
                int new_x = here.X + model->items.at(i_i).second[i_v].first;
                int new_y = here.Y + model->items.at(i_i).second[i_v].second;
                if( ( new_x >= here.d1 ) || ( new_x < 0 ) || ( new_y >= here.d2 ) || ( new_y < 0 ) )
                    continue;
                count++;
                GridSite there = GridSite( new_x, new_y, *grid );

                // GridSite there = GridSite( i_g, *grid )
                //     + GridSite(
                //         model->items.at(i_i).second[i_v].first,
                //         model->items.at(i_i).second[i_v].second,
                //         *grid
                //         );
                gshere  = ( grid->operator[](  here ) == GridSite::Defect ? 1 : -1 );
                gsthere = ( grid->operator[]( there ) == GridSite::Defect ? 1 : -1 );
                values->at( i_i ) += gshere * gsthere;
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

Expospaced::Expospaced( double base, int maxval ) {
    int current_exp = 0;
    int current_d = 1;
    do {

        items.push_back( make_pair( current_d, vector< pair<int, int> >()  ) );
        for( int v1 = - current_d - 1; v1 <= current_d + 1; v1 ++ ) {
            for( int v2 = - current_d - 1; v2 <= current_d + 1; v2 ++ ) {
                if( int( round( sqrt( v1 * v1 + v2 * v2 ) ) ) == current_d ) {
                    items.back().second.push_back( make_pair( v1, v2 ) );
                }
            }
        }

        while( int( round( pow( base, current_exp ) ) ) <= current_d )
            current_exp++;
        current_d = int( round( pow( base, current_exp ) ) );
    } while( current_d < maxval );

    keyname = "Expospaced";
}

OrtExpospaced::OrtExpospaced( double base, int maxval ) {
    int current_exp = 0;
    int current_d = 1;
    do {

        items.push_back( make_pair( current_d, vector< pair<int, int> >()  ) );
        items.back().second.push_back( make_pair( current_d, 0 ) );
        // items.back().second.push_back( make_pair(-current_d, 0 ) );
        items.back().second.push_back( make_pair( 0, current_d ) );
        // items.back().second.push_back( make_pair( 0,-current_d ) );

        while( int( round( pow( base, current_exp ) ) ) <= current_d )
            current_exp++;
        current_d = int( round( pow( base, current_exp ) ) );
    } while( current_d < maxval );

    keyname = "OrtExpospaced";
}