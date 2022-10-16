#include "Replicator.h"

using namespace std;

string ReplicatorParams::to_string() {
    stringstream ss;
    string sep = " | ";
    string trimmed_savepath = (
        save_path.length() < 13 ?
        string( 13 - save_path.length(), ' ' ) + save_path :
        save_path.substr( save_path.length() - 13 ) );

    ss << date::format("%H:%M", chrono::system_clock::now()) << sep
       << trimmed_savepath << sep
       << setw( 4 ) << side << sep 
       << setw( 12 )<< defects_frac << sep
       << setw( 5 ) << gamma << sep 
       << setw( 7 ) << chunk_size << sep
       << setw( 9 ) << tolerance << sep
       << (CF_model ? CF_model->keyname : "Disabled  ") << sep
       << (to_deposit ? to_deposit->keyname : "Disabled") << sep
       << (percolation ? "Both" : "Disabled") << sep;

    return ss.str();
}

string ReplicatorParams::header() {
    return "TIME  | PATHNAME      | SIDE | DEFECTS_FRAC | GAMMA | N_CHUNK | TOLERANCE | CORR_RANGE | POLYMERS | PERCOLATION";
}

double ReplicatorParams::size() const {
    return side*side;
}

ReplicatorThread::ReplicatorThread(Replicator* thrower, int id)
    : thrower( thrower ), id( id ),
    g( thrower->params.side ), fcg( thrower->params.side, thrower->params.gamma ),
    gfp( thrower->params.to_deposit == nullptr ? 0 : thrower->params.to_deposit->N, thrower->params.size() )  {

    CFcalc  = nullptr;
    CFcalcH = nullptr;

    perc = nullptr;
    if( thrower->params.percolation )
        perc = new Percolator( g );
}

ReplicatorThread::~ReplicatorThread() {
    delete CFcalc;
    delete CFcalcH;
    delete perc;
}

/*static*/ void ReplicatorThread::thread_worker (ReplicatorThread* data) {
    data->thrower->mux->lock();
    while( data->thrower->runned_replicas < data->thrower->total_replicas_to_run ) {
        int current_replica = ++ data->thrower->runned_replicas;
        if( data->thrower->params.verbose )
            if( current_replica % max( data->thrower->total_replicas_to_run / 100, 1 ) == 0 )
                cout<< "[prg] " << current_replica << '/' << data->thrower->total_replicas_to_run << "\t" << " [" << data->id << "]                          \r" <<flush;
        data->thrower->mux->unlock();

        // Clean the grid
        GridFiller::clean( data->g );

        // Populate h field
        data->fcg.repopulate();

        // Deposit defects
        GridFiller::ranked_insertion(
            data->g,
            data->fcg.h,
            data->thrower->params.defects_frac * data->thrower->params.side * data->thrower->params.side
            );

        // Compute correlation functions
        if( data->thrower->params.CF_model ) {
            data->thrower->update_NewCF_averages(
                data->CFcalc->calculate()
            );
            data->thrower->update_NewCFH_averages(
                data->CFcalcH->calculate()
            );
        }

        // Deposit polymers
        if( data->thrower->params.to_deposit != nullptr ) {
            double occupied_sites = data->gfp.fill( data->g, *data->thrower->params.to_deposit );
            data->thrower->update_dep_averages( occupied_sites );
        }

        // Compute percolation infos
        if( data->thrower->params.percolation ) {
            data->thrower->update_perc_averages(
                data->perc->is_percolating( GridSite::Defect ),
                data->perc->is_percolating( GridSite::Atom )
            );
        }

        // Print the grid
        if( data->thrower->params.draw ) {
            data->g.print_data( data->thrower->params.save_path + "/draw_" + to_string(current_replica) + ".txt" );
        }

        // Proceed to other replicas
        data->thrower->mux->lock();
    }
    data->thrower->mux->unlock();
}

void ReplicatorThread::start() {
    thrd = new thread( ReplicatorThread::thread_worker, this );
}

void ReplicatorThread::join() {
    thrd->join();
    delete thrd;
}

void Replicator::addChunk() {
    total_replicas_to_run *= 2;
    if( total_replicas_to_run == 0 )
        total_replicas_to_run = params.chunk_size;
    
    fills->reserve( total_replicas_to_run );
    
    for( ReplicatorThread* rt : ongoing )
        rt->start();
    for( ReplicatorThread* rt : ongoing )
        rt->join();
}

Replicator::Replicator( ReplicatorParams suggested_params ) :
    params( suggested_params ), mux( new mutex() ) {

    // Compute save path
    // - Ensure uniqueness adding a progressive
    int progressive = 0;
    while( filesystem::exists( params.save_path + "/f_" + to_string( progressive ) ) )
        progressive++;
    params.save_path = params.save_path + "/f_" + to_string( progressive );
    filesystem::create_directory( params.save_path );

    // Create correlation function structures
    CF_avg = nullptr;
    CFH_avg = nullptr;
    if( params.CF_model ) {
        CF_avg = new vector< double > ( params.CF_model->items.size(), 0 );
        CFH_avg = new vector< double > ( params.CF_model->items.size(), 0 );
    }

    // Initialize occupation stats
    fills = new vector< unsigned int >();

    // Initialize percolation stats
    defperc_count = 0;
    atmperc_count = 0;

    // Initialize the thread managers
    for( int i=0; i < params.n_threads; i++ ) {
        ongoing.push_back( new ReplicatorThread( this, i ) );
    }

    // Initialize the replicas counter
    runned_replicas = 0;
    total_replicas_to_run = 0;    
}

Replicator::~Replicator() {
    delete mux;
    for( ReplicatorThread* r : ongoing )
        delete r;
    delete CF_avg;
}

void Replicator::update_NewCF_averages( const vector< double >* cf ) {
    mux->lock();
    for( int i = 0; i < cf->size(); i++ ) {
        CF_avg->at( i ) += cf->at( i );
    }
    mux->unlock();
}

void Replicator::update_NewCFH_averages( const vector< double >* cf ) {
    mux->lock();
    for( int i = 0; i < cf->size(); i++ ) {
        CFH_avg->at( i ) += cf->at( i );
    }
    mux->unlock();
}

void Replicator::update_dep_averages( double fill_frac ) {
    mux->lock();
    if( fill_frac > UINT_MAX ) {
        cout<<"Warning! Overflow exeption!";
    }
    fills->push_back( fill_frac );
    mux->unlock();
}

double Replicator::fill_avg( unsigned int threshold ) const {
    if( runned_replicas < threshold )
        threshold = runned_replicas;
    double result = 0;
    for( int i = 0; i < threshold; i++ )
        result += fills->at(i);
    return result / threshold;
}

double Replicator::fill_std_fromnp( unsigned int threshold ) const {
    if( runned_replicas < threshold )
        threshold = runned_replicas;

    double avg = fill_avg( threshold );
    double result = 0;

    for( int i = 0; i < threshold; i++ )
        result += ( ( fills->at(i) - avg ) * ( fills->at(i) - avg ) );

    return sqrt( result / ( threshold - 1 ) );
}

double Replicator::fill_std_fromgauss( unsigned int threshold ) const {
    if( runned_replicas < threshold )
        threshold = runned_replicas;
    
    vector< int > counts( params.size(), 0 );
    for( int i = 0; i < counts.size(); i++ ) {
        counts[i] = 0;
    }
    for( int i = 0; i < threshold; i++ ) {
        counts[ fills->at(i) ]++;
    }

    vector< double > x;
    vector< double > y;
    for( int i = 0; i < counts.size(); i++ ) {
        if( counts[i] > 0 ) {
            x.push_back( i );
            y.push_back( counts[i] );
        }
    }

    if( x.size() < 3 ) return 0;

    double avg = fill_avg( threshold );
    auto params = curve_fit( gaussian, { (double)threshold, avg, avg * 0.1 }, x, y );

    return abs( params[2] );
}

double Replicator::fill_std( unsigned int threshold ) const {
    if( WHICHSTD == STD_FROMFORMULA )
        return fill_std_fromnp( threshold );
    else if( WHICHSTD == STD_FROMGAUSS )
        return fill_std_fromgauss( threshold );
    return 0;
}

void Replicator::update_perc_averages( bool def_perc, bool atm_perc ) {
    mux->lock();
    defperc_count += ( def_perc ? 1 : 0 );
    atmperc_count += ( atm_perc ? 1 : 0 );
    mux->unlock();
}

void Replicator::run() {

    ofstream out_chunks( params.save_path + "/chunks.txt");
    out_chunks << " RunnedReplicas | Std | Variation | Items"<<endl;
    if( params.verbose )
        cout<< "RunnedReplicas | Std | Variation " <<endl;

    vector< double > stds;
    double variation = INT_MAX;
    while(1) {
        addChunk();

        if( params.to_deposit == nullptr ) {
            break;
        }

        if( params.tolerance < 0 ) {
            break;
        }

        if( runned_replicas > MAX_REPS ) {
            break;
        }

        stds.push_back( fill_std() );

        if( stds.size() > 2 ) {
            variation = *max_element( stds.end() - 2, stds.end() ) / *min_element( stds.end() - 2, stds.end() );
        }

        out_chunks << runned_replicas << ", " << stds.back() << ", " << variation << ", [";
        for( unsigned int d : *fills )
            out_chunks<<d<<',';
        out_chunks<<"]"<<endl;
        if( params.verbose )
            cout << runned_replicas << ", " << stds.back() << " (" << fill_std() << "), " << variation << endl;

        if( isnan( variation ) )
            break;
        if( variation < params.tolerance + 1 )
            break;
    }
}

void Replicator::save_data() {
    // Print details
    ofstream out_details( params.save_path + "/details.txt");
    out_details<<"{\n\"side\":\t"<<params.side
               <<",\n\"defects_frac\":\t"<<params.defects_frac
               <<",\n\"gamma\":\t"<<params.gamma
               <<",\n\"chunk_size\":\t"<<params.chunk_size
               <<",\n\"tolerance\":\t"<<params.tolerance
               <<",\n\"runned_replicas\":\t"<<runned_replicas
               <<",\n\"CF_Model\":\t\""<<(params.CF_model ? params.CF_model->keyname : "") << "\""
               <<",\n\"dep_polymers\":\t\""<<(params.to_deposit ? params.to_deposit->keyname : "") << "\""
               <<",\n\"save_path\":\t\"" << params.save_path << "\""
               <<",\n\"n_threads\":\t\"" << params.n_threads << "\""
               <<",\n\"draw\":\t" << (params.draw ? "True" : "False" ) << ""
               <<"\n}"<<endl;
    out_details.close();

    // Print correlation function stuff
    if( params.CF_model ) {
        ofstream out_corr( params.save_path + "/CF_avg.txt");
        for( int i=0; i < CF_avg->size(); i++ ) {
            out_corr << params.CF_model->items[i].first << '\t' << CF_avg->at( i ) / runned_replicas << '\t' << params.CF_model->items[i].second.size() << '\n';
        }
        out_corr.close();

        ofstream out_corr_H( params.save_path + "/CFH_avg.txt");
        for( int i=0; i < CFH_avg->size(); i++ ) {
            out_corr_H << params.CF_model->items[i].first << '\t' << CFH_avg->at( i ) / runned_replicas << '\t' << params.CF_model->items[i].second.size() << '\n';
        }
        out_corr_H.close();
    }

    // Print polymers deposition results
    if( params.to_deposit ) {
        double avg = fill_avg();
        double std = fill_std();

        ofstream out_deposition( params.save_path + "/deposition.txt");
        out_deposition<<"{\n\"dep_polymers\":\t\""<<params.to_deposit->keyname<<"\""
                   <<",\n\"occupation_average\":\t"<< avg
                   <<",\n\"occupation_std\":\t"<<std
                   <<",\n\"occupation_fraction_average\":\t"<< avg / params.size()
                   <<",\n\"occupation_fraction_std\":\t"<<std / params.size()
                   <<",\n\"pj_over_1_minus_q_avg\":\t"<< avg / params.size() / ( 1 - params.defects_frac )
                   <<",\n\"pj_over_1_minus_q_std\":\t"<< std / params.size() / ( 1 - params.defects_frac )
                   <<",\n\"occupation_history\":\t[";
        for( unsigned int d : *fills )
            out_deposition<<d<<',';
        out_deposition<<"]\n}"<<endl;
        out_deposition.close();
    }

    // Print percolation results
    if( params.percolation ) {
        double defperc_avg = defperc_count / runned_replicas;
        double atmperc_avg = atmperc_count / runned_replicas;

        ofstream out_deposition( params.save_path + "/percolation.txt");
        out_deposition<<"{"
                   <<" \n\"defperc_count\":\t"<<defperc_count
                   <<",\n\"defperc_avg\":\t"<<defperc_avg
                   <<",\n\"atmperc_count\":\t"<<atmperc_count
                   <<",\n\"atmperc_avg\":\t"<<atmperc_avg
                   <<"\n}"<<endl;
        out_deposition.close();
    }
}

int32_t Replicator::total_runned() {
    return runned_replicas;
}