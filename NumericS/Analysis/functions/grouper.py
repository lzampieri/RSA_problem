import numpy as np

def extract_keys( data, key_list ):
    return { k:v for k,v in data.items() if k in key_list }

def hash_dict( data ):
    return '_'.join( [ str( v ) for v in data.values() ] )

def extract_xy( data, keys, x, y ):
    output = {}
    for d in data:
        values = extract_keys( d, keys )
        key = hash_dict( values )
        output.setdefault( key, {} ).update( values )
        output[key].setdefault( 'x', [] ).append( x( d ) )
        output[key].setdefault( 'y', [] ).append( y( d ) )
    
    return np.array( list( output.values() ) )
