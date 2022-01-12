#include "Replicator.h"

using namespace std;

Replicator::Replicator( int side, double defects_frac, double gamma, int n_replies, std::string save_path ) :
    side( side ), defects_frac( defects_frac ), gamma( gamma ), n_replies( n_replies ),
    corr_range( -1 ), fcg( side ), g( side ) {
    
    if( save_path == "" )
        this->save_path = date::format("%Y%m%d", chrono::system_clock::now());
    else
        this->save_path = save_path;
}

Replicator::~Replicator() {
    delete CF_H;
    delete CF_D;
}

void Replicator::enable_correlators( int corr_range ) {
    this->corr_range = min( corr_range, side / 2 );

    CF_H = new CorrFunc::Expospaced<double>( &fcg.h, 1.4, this->corr_range, 1 );
    CF_D = new CorrFunc::Expospaced<int>   ( &g    , 1.4, this->corr_range, 1 );
}

void Replicator::run_replica( vector< double >* CF_H_avg, vector< double >* CF_D_avg ) {

    GridFiller::clean(g);

    GridFiller::iid( fcg.h );
    fcg.fourier_transform();
    fcg.multiply_fft_new( gamma );
    fcg.fourier_transform(true);

    GridFiller::ranked_insertion( g, fcg.h, defects_frac * side * side );
    
    if( corr_range > 0 ) { // If correlation calcolation is enabled
        const vector< CorrFunc::Datapoint >* cfh = CF_H->compute_corr_function ( );
        const vector< CorrFunc::Datapoint >* cfd = CF_D->compute_corr_function ( );

        for( int i=0; i < CF_H_avg->size(); i++ ) {
            CF_H_avg->at(i) += cfh->at(i).value;
        }
        
        for( int i=0; i < CF_D_avg->size(); i++ ) {
            CF_D_avg->at(i) += cfd->at(i).value;
        }
    }

}

string Replicator::run() {
    // String for indexing
    stringstream ss;
    string sep = " | ";

    // Compute save path
    int progressive = 0;
    while( filesystem::exists( save_path + "_" + to_string( progressive ) ) )
        progressive++;
    string actual_path = save_path + "_" + to_string( progressive );
    if( !filesystem::exists( actual_path ) )
        filesystem::create_directory( actual_path );
    ss << actual_path << '\t' << sep;

    // Print details
    ofstream out_details( actual_path + "/details.txt");
    out_details<<"{\n\"side\":\t"<<side<<",\n\"defects_frac\":\t"<<defects_frac<<
               ",\n\"gamma\":\t"<<gamma<<",\n\"replies\":"<<n_replies<<
               ",\n\"corr_range\":\t"<<corr_range<<"\n}"<<endl;
    ss << setw( 4 ) << side << sep 
       << setw( 12 )<< defects_frac << sep
       << setw( 5 ) << gamma << sep 
       << setw( 9 ) << n_replies << sep 
       << setw( 10 )<< corr_range << sep;

    // Correlation function management
    vector< double >* CF_H_avg = nullptr;
    vector< double >* CF_D_avg = nullptr;
    if( corr_range > 0 ) {// If correlation calcolation is enabled
        CF_H_avg = new vector< double > ( CF_H->is->size(), 0 );
        CF_D_avg = new vector< double > ( CF_D->is->size(), 0 );
    }

    // Run replicas
    cout<<"[log] Starting replicas performing..."<<endl;
    for( int i=0; i < n_replies; i++ ) {
        run_replica( CF_H_avg, CF_D_avg );
        cout<< "[prg] " << i << "\t/" << n_replies << '\r' <<flush;
    }
    cout<<"[log] End replicas.  "<<endl;

    // Print correlation function stuff
    if( corr_range > 0 ) {
        ofstream out_corr( actual_path + "/CF_H_avg.txt");
        for( int i=0; i < CF_H->is->size(); i++ ) {
            CF_H_avg->at(i) /= n_replies;
            out_corr<< CF_H->is->at(i) << '\t' << CF_H_avg->at(i) << '\n';
        }

        ofstream out2_corr( actual_path + "/CF_D_avg.txt");
        for( int i=0; i < CF_D->is->size(); i++ ) {
            CF_D_avg->at(i) /= n_replies;
            out_corr<< CF_D->is->at(i) << '\t' << CF_D_avg->at(i) << '\n';
        }
    }

    // Clean
    delete CF_H_avg;
    delete CF_D_avg;

    // Return a string in the form PATHNAME | SIDE | DEFECTS_FRAC | GAMMA | N_REPLIES | CORR_RANGE
    return ss.str();

}

std::string Replicator::run_replicator( ReplicatorParams params ) {
    Replicator r( params.side, params.defects_frac, params.gamma, params.n_replies );
    if( params.corr_range > 0 ) r.enable_correlators( params.corr_range );
    return r.run();
}