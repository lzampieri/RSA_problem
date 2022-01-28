#include "Replicator.h"

using namespace std;

std::string ReplicatorParams::to_string() {
    stringstream ss;
    string sep = " | ";
    ss << save_path << '\t' << sep;
    ss << setw( 4 ) << side << sep 
       << setw( 12 )<< defects_frac << sep
       << setw( 5 ) << gamma << sep 
       << setw( 9 ) << n_replies << sep 
       << (CF_model ? CF_model->keyname : "") << sep
       << (to_deposit ? to_deposit->keyname : "") << sep
       << (percolation ? "Percolation" : "") << sep;;

    return ss.str();
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
    while( data->thrower->replicas_to_run > 0 ) {
        int current_replica = data->thrower->replicas_to_run;
        data->thrower->replicas_to_run --;
        if( current_replica % max( data->thrower->params.n_replies / 10, 1 ) == 0 )
            cout<< "[prg] " << ( data->thrower->params.n_replies - current_replica ) * 100 / data->thrower->params.n_replies << "%\t" << " [" << data->id << "]                          \r" <<flush;
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
            int occupied_sites = GridFiller::fillWithPolymers( data->g, *data->thrower->params.to_deposit );
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

Replicator::Replicator( ReplicatorParams suggested_params ) :
    params( suggested_params ), mux( new mutex() ) {

    // Compute save path
    // - If no other requirements, define from date
    if( params.save_path == "" )
        params.save_path = date::format("%Y%m%d", chrono::system_clock::now());
    // - Ensure uniqueness adding a progressive
    int progressive = 0;
    while( filesystem::exists( params.save_path + "_" + to_string( progressive ) ) )
        progressive++;
    params.save_path = params.save_path + "_" + to_string( progressive );
    if( !filesystem::exists( params.save_path ) )
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
    replicas_to_run = params.n_replies;
    
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

void Replicator::update_dep_averages( int occupied_sites ) {
    mux->lock();
    fillfrac_sum  += occupied_sites;
    fillfrac_sum2 += occupied_sites*occupied_sites;
    mux->unlock();
}

void Replicator::update_perc_averages( bool def_perc, bool atm_perc ) {
    mux->lock();
    defperc_count += ( def_perc ? 1 : 0 );
    atmperc_count += ( atm_perc ? 1 : 0 );
    mux->unlock();
}

void Replicator::run() {
    cout<<"[log] Start: " << params.to_string() <<endl;
    for( ReplicatorThread* rt : ongoing )
        rt->start();
    for( ReplicatorThread* rt : ongoing )
        rt->join();
}

string Replicator::save_data() {
    // Print details
    ofstream out_details( params.save_path + "/details.txt");
    out_details<<"{\n\"side\":\t"<<params.side
               <<",\n\"defects_frac\":\t"<<params.defects_frac
               <<",\n\"gamma\":\t"<<params.gamma
               <<",\n\"replies\":\t"<<params.n_replies
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
            out_corrH<< params.CF_model->is[i] << '\t' << CF_H_avg->at(i) / params.n_replies << '\n';
            out_corrD<< params.CF_model->is[i] << '\t' << CF_D_avg->at(i) / params.n_replies << '\n';
        }
        out_corrH.close();
        out_corrD.close();
    }

    // Print polymers deposition results
    if( params.to_deposit ) {
        double avg = fillfrac_sum / params.n_replies;
        double std = sqrt( ( fillfrac_sum2 - fillfrac_sum*fillfrac_sum/params.n_replies ) / ( params.n_replies - 1 ) );

        ofstream out_deposition( params.save_path + "/deposition.txt");
        out_deposition<<"{\n\"dep_polymers\":\t\""<<params.to_deposit->keyname<<"\""
                   <<",\n\"occupation_average\":\t"<<avg
                   <<",\n\"occupation_std\":\t"<<std
                   <<",\n\"occupation_fraction_average\":\t"<< avg / (params.side*params.side)
                   <<",\n\"occupation_fraction_std\":\t"<<std / (params.side*params.side)
                   <<"\n}"<<endl;
        out_deposition.close();
    }

    // Print percolation results
    if( params.percolation ) {
        double defperc_avg = defperc_count / params.n_replies;
        double atmperc_avg = atmperc_count / params.n_replies;

        ofstream out_deposition( params.save_path + "/percolation.txt");
        out_deposition<<"{"
                   <<" \n\"defperc_count\":\t"<<defperc_count
                   <<",\n\"defperc_avg\":\t"<<defperc_avg
                   <<",\n\"atmperc_count\":\t"<<atmperc_count
                   <<",\n\"atmperc_avg\":\t"<<atmperc_avg
                   <<"\n}"<<endl;
        out_deposition.close();
        cout<<"Percolation averages:\nDefects: "<<defperc_avg<<"\nAtoms: "<<atmperc_avg<<endl; //todo remove
    }

    return params.to_string();
}