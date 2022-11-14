#include "AutoScanner.h"

using namespace std;

string AutoScanner::populate() {
    string folder = computeFolder();

    for( int s : sides )
        for( double g : gammas )
            for( double q : qs )
                for( PolymersFactory* p : ps ) {
                    rps.push_back( ReplicatorParams(
                     // Size DefectsFracs Gamma ChunkSize   Tolerance  CFModel  Polymers                                   NThreads,  Draw  Verbose  SavePath
                        s,   q,           g,    chunk_size, tolerance, CFmodel, ( p == nullptr ? nullptr : p->create(s) ), n_threads, draw, verbose, folder
                    ));
                }

    return folder;
}

AutoScanner::AutoScanner() {

    // Default params
    chunk_size = 10000;
    tolerance = -1;

    sides = { 512 };
    gammas = { 0.4, 0.8, 1.2, 1.6 };
    qs = { 0.5 };
    
    for( PolymersFactory* p : PolymersFactory::StdPolymers ) {
        ps.push_back( p );
    }

    draw = false;
    verbose = true;
    n_threads = 32;

    CFmodel = nullptr;

}

void AutoScanner::loadFromTxt( string file ) {
    
    filename = file;
    
    ifstream in( file );
    string temp, content;

    while( smart_getline( in, temp ) ) {
        
        content = temp.substr( 2 );

        switch( temp.at( 0 ) ){
            case 's':
                chunk_size = stoi( content );
                break;

            case 't':
                tolerance = stod( content );
                break;

            case 'L':
                string_to_array( content, sides );
                break;

            case 'g':
                string_to_array( content, gammas );
                break;

            case 'q':
                string_to_array( content, qs );
                break;

            case 'p':
                string_to_array( content, ps );
                break;

            case 'd':
                draw = stoi( content );
                break;

            case 'v':
                verbose = stoi( content );
                break;

            case 'T':
                n_threads = stoi( content );
                break;

            default:
                throw invalid_argument("Invalid sintax.");

        }
    }
}

string AutoScanner::computeFolder( ) {

    string folder = filename;

    // Ensure only numbers and letters
    folder.erase( remove_if( folder.begin(), folder.end(),
        [](char c) { return !isalpha(c) && !isdigit(c); } ),
        folder.end()
        );

    // Append current date and versioning
    folder += "_" + date::format("%Y%m%d", chrono::system_clock::now()) + "_v2";

    // Check for existence
    if( filesystem::exists( folder ) ) {
        int progressive = 0;
        // If already exists, append progressive
        while( filesystem::exists( folder + "_" + to_string( progressive ) ) )
            progressive ++;
        folder += "_" + to_string( progressive );
    }

    // Create
    filesystem::create_directory( folder );

    // Save there a copy of the configuration file
    saveToTxt( folder + "/params.scan" );

    return folder;
}

void AutoScanner::saveToTxt( string file ) {
    ofstream out( file );

    out << "s " << chunk_size << '\n'
        << "t " << tolerance << '\n'
        << "L " << array_to_string( sides ) << '\n'
        << "g " << array_to_string( gammas ) << '\n'
        << "q " << array_to_string( qs ) << '\n'
        << "p " << array_to_string( ps ) << '\n'
        << "d " << draw << '\n'
        << "v " << verbose << '\n'
        << "T " << n_threads << endl;

    out.close();
}

void AutoScanner::string_to_array( const string text, vector<double>& vec ) {
    stringstream ss( text );
    vec.clear();
    double t;
    while( ss >> t )
        vec.push_back( t );
}

void AutoScanner::string_to_array( const string text, vector<int>& vec ) {
    stringstream ss( text );
    vec.clear();
    int t;
    while( ss >> t )
        vec.push_back( t );
}

void AutoScanner::string_to_array( const string text, vector<PolymersFactory*>& vec ) {
    stringstream ss( text );
    string temp;
    vec.clear();
    while( ss >> temp ) {
        for( PolymersFactory* p : PolymersFactory::StdPolymers ) {
            if( p->factname() == temp ) {
                vec.push_back( p );
                break;
            }
        }
    }
}

string AutoScanner::array_to_string( const std::vector<double>& vec ) {
    stringstream ss( "" );
    for( double t : vec )
        ss << t << ' ';
    return ss.str();
}

string AutoScanner::array_to_string( const std::vector<int>& vec ) {
    stringstream ss( "" );
    for( int t : vec )
        ss << t << ' ';
    return ss.str();
}

string AutoScanner::array_to_string( const std::vector<PolymersFactory*>& vec ) {
    stringstream ss( "" );
    for( const PolymersFactory* t : vec )
        ss << ( t == nullptr ? "Null ptr" : t->factname() ) << ' ';
    return ss.str();
}

bool AutoScanner::smart_getline( ifstream& in, string& str ) {
    do {

        if( in.eof() ) return false;
        getline( in, str );

    } while( str.length() == 0 );

    return true;
}