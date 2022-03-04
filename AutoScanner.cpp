#include "AutoScanner.h"

using namespace std;

void AutoScanner::populate() {
    for( int s : sides )
        for( double g : gammas )
            for( double q : qs )
                for( PolymersFactory* p : ps )
                    rps.push_back( ReplicatorParams(
                     // Size DefectsFracs Gamma ChunkSize   Tolerance  CFModel  Polymers      Percolation  NThreads,  Draw  Verbose SavePath
                        s,   q,           g,    chunk_size, tolerance, CFmodel, p->create(s), percolation, n_threads, draw, verbose
                    ));
}

void AutoScanner::loadFromTxt( string file ) {
    ifstream in( file );
    string temp;

    in  >> chunk_size
        >> tolerance;
    
    smart_getline( in, temp );
    string_to_array( temp, sides );
    
    smart_getline( in, temp );
    string_to_array( temp, gammas );

    smart_getline( in, temp );
    string_to_array( temp, qs );

    smart_getline( in, temp );
    string_to_array( temp, ps );
    
    in >> percolation
        >> draw
        >> verbose
        >> n_threads;
}

void AutoScanner::saveToTxt( string file ) {
    ofstream out( file );

    out << chunk_size << '\n'
        << tolerance << '\n'
        << array_to_string( sides ) << '\n'
        << array_to_string( gammas ) << '\n'
        << array_to_string( qs ) << '\n'
        << array_to_string( ps ) << '\n'
        << percolation << '\n'
        << draw << '\n'
        << verbose << '\n'
        << n_threads << endl;

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
    while( ss >> temp ) {
        for( PolymersFactory* p : PolymersFactory::StdPolymers ) {
            if( p->factname() == temp ) {
                ps.push_back( p );
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
        ss << t->factname() << ' ';
    return ss.str();
}

void AutoScanner::smart_getline( ifstream& in, string& str ) {
    do {
        getline( in, str );
    } while( str.length() == 0 );
}