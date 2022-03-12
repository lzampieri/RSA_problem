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
    g( thrower->params.side ), fcg( thrower->params.side )  {

    CF_H = nullptr;
    CF_D = nullptr;
    if( thrower->params.CF_model ) {
        CF_H = new CorrFunc::Calculator<double>( thrower->params.CF_model, &fcg.h );
        CF_D = new CorrFunc::Calculator<int>   ( thrower->params.CF_model, &g     );
    }

    perc = nullptr;
    if( thrower->params.percolation )
        perc = new Percolator( g );
}

ReplicatorThread::~ReplicatorThread() {
    delete CF_H;
    delete CF_D;
    delete perc;
}

/*static*/ void ReplicatorThread::thread_worker (ReplicatorThread* data) {
    data->thrower->mux->lock();
    while( data->thrower->runned_replicas < data->thrower->total_replicas_to_run ) {
        int current_replica = ++ data->thrower->runned_replicas;
        if( data->thrower->params.verbose )
            if( current_replica % max( data->thrower->params.chunk_size / 100, 1 ) == 0 )
                cout<< "[prg] " << current_replica << '/' << data->thrower->params.chunk_size << "\t" << " [" << data->id << "]                          \r" <<flush;
        data->thrower->mux->unlock();

        // Clean the grid
        GridFiller::clean( data->g );

        // Populate defects
        GridFiller::iid( data->fcg.h );
        data->fcg.fourier_transform();
        data->fcg.multiply_fft_new( data->thrower->params.gamma );
        data->fcg.fourier_transform(true);
        GridFiller::ranked_insertion(
            data->g,
            data->fcg.h,
            data->thrower->params.defects_frac * data->thrower->params.side * data->thrower->params.side
            );

        // Compute correlation functions
        if( data->thrower->params.CF_model ) {
            const vector< CorrFunc::Datapoint >* cfh = data->CF_H->compute_corr_function ( );
            const vector< CorrFunc::Datapoint >* cfd = data->CF_D->compute_corr_function ( );
            data->thrower->update_CF_averages( cfh, cfd );
            
        }

        // Deposit polymers
        if( data->thrower->params.to_deposit != nullptr ) {
            double occupied_sites = GridFiller::fillWithPolymers( data->g, *data->thrower->params.to_deposit );
            data->thrower->update_dep_averages( occupied_sites / data->thrower->params.size() );
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
    total_replicas_to_run += params.chunk_size;
    
    for( ReplicatorThread* rt : ongoing )
        rt->start();
    for( ReplicatorThread* rt : ongoing )
        rt->join();

    fillfrac_stds.push_back( fillfrac_std() );
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
    CF_H_avg = nullptr;
    CF_D_avg = nullptr;
    if( params.CF_model ) {
        CF_D_avg = new vector< double > ( params.CF_model->is.size(), 0 );
        CF_H_avg = new vector< double > ( params.CF_model->is.size(), 0 );
    }

    // Initialize occupation stats
    fillfrac_sum = 0;
    fillfrac_sum2= 0;

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
    delete CF_D_avg;
    delete CF_H_avg;
}

void Replicator::update_CF_averages( const std::vector< CorrFunc::Datapoint >* cfh, const std::vector< CorrFunc::Datapoint >* cfd ) {
    mux->lock();
    for( int i=0; i < CF_H_avg->size(); i++ ) {
        CF_H_avg->at(i) += cfh->at(i).value;
    }
    for( int i=0; i < CF_D_avg->size(); i++ ) {
        CF_D_avg->at(i) += cfd->at(i).value;
    }
    mux->unlock();
}

void Replicator::update_dep_averages( double fill_frac ) {
    mux->lock();
    fillfrac_sum  += fill_frac;
    fillfrac_sum2 += fill_frac*fill_frac;
    mux->unlock();
}

double Replicator::fillfrac_std() const {
    return sqrt( fillfrac_sum2 / runned_replicas - pow(fillfrac_sum / runned_replicas, 2) );
}

void Replicator::update_perc_averages( bool def_perc, bool atm_perc ) {
    mux->lock();
    defperc_count += ( def_perc ? 1 : 0 );
    atmperc_count += ( atm_perc ? 1 : 0 );
    mux->unlock();
}

void Replicator::run() {
    while(1) {
        addChunk();
        if( params.to_deposit && fillfrac_stds.size() > 1 ) {
            double last = fillfrac_stds[ fillfrac_stds.size() - 1 ];
            double prelast = fillfrac_stds[ fillfrac_stds.size() - 2 ];
            if( params.verbose )
                cout << "Variation: " << abs( last - prelast ) / min( last, prelast ) << endl;
            if( abs( last - prelast ) / min( last, prelast ) < params.tolerance )
                break;
        }
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
               <<",\n\"runned_chunks\":\t"<< (runned_replicas+1) / params.chunk_size
               <<",\n\"CF_Model\":\t\""<<(params.CF_model ? params.CF_model->keyname : "") << "\""
               <<",\n\"dep_polymers\":\t\""<<(params.to_deposit ? params.to_deposit->keyname : "") << "\""
               <<",\n\"save_path\":\t\"" << params.save_path << "\""
               <<",\n\"n_threads\":\t\"" << params.n_threads << "\""
               <<",\n\"draw\":\t" << (params.draw ? "True" : "False" ) << ""
               <<"\n}"<<endl;
    out_details.close();

    // Print correlation function stuff
    if( params.CF_model ) {
        ofstream out_corrH( params.save_path + "/CF_H_avg.txt");
        ofstream out_corrD( params.save_path + "/CF_D_avg.txt");
        for( int i=0; i < params.CF_model->is.size(); i++ ) {
            out_corrH<< params.CF_model->is[i] << '\t' << CF_H_avg->at(i) / runned_replicas << '\n';
            out_corrD<< params.CF_model->is[i] << '\t' << CF_D_avg->at(i) / runned_replicas << '\n';
        }
        out_corrH.close();
        out_corrD.close();
    }

    // Print polymers deposition results
    if( params.to_deposit ) {
        double avg = fillfrac_sum / runned_replicas;
        double std = fillfrac_std();

        ofstream out_deposition( params.save_path + "/deposition.txt");
        out_deposition<<"{\n\"dep_polymers\":\t\""<<params.to_deposit->keyname<<"\""
                   <<",\n\"occupation_fraction_average\":\t"<< avg
                   <<",\n\"occupation_fraction_std\":\t"<<std
                   <<",\n\"pj_over_1_minus_q_avg\":\t"<< avg / ( 1 - params.defects_frac )
                   <<",\n\"pj_over_1_minus_q_std\":\t"<< std / ( 1 - params.defects_frac )
                   <<",\n\"occupation_fraction_std_history\":\t[";
        for( double d : fillfrac_stds )
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

int Replicator::runned_chunks() {
    return ( runned_replicas + 1 ) / params.chunk_size;
}