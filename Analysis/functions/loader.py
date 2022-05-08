import numpy as np
import json
from pathlib import Path
import os
from glob import glob
import re

# Load single file
def load_file( filename ):
    text = filename.read_text()
    text = text.replace("nan","0")
    text = text.replace("False", "false" )
    text = text.replace("True", "true" )
    text = text.replace(",]", "]" )
    return json.loads( text )

def load_data( ):
    data = []

    for file in glob("../*Analysis/**/details.txt", recursive=True):
        item = {}
        d = Path( os.path.dirname( file ) )

        # Verify that everything exists
        if( not (
            ( d / 'details.txt').exists() and
            ( d / 'deposition.txt').exists() and
            ( d / 'percolation.txt').exists() and
            ( d / 'chunks.txt').exists() ) ):
            continue

        # Load data
        item.update( load_file( d / 'details.txt'     ) )
        item.update( load_file( d / 'deposition.txt'  ) )
        item.update( load_file( d / 'percolation.txt' ) )
        with open( d / "chunks.txt" ) as file:
            lines = file.readlines()
            chunks = json.loads( "[" + lines[-1].replace(",]","]").replace("nan","0") + "]" )[-1]
            item['chunks'] = chunks


        # Prepare histogram
        bins = np.arange( min( chunks ), max( chunks ) + 1 ) - 0.5
        h_all, _ = np.histogram( chunks, bins= bins, density=True )
        h_x_all = ( bins[:-1] + 0.5 ) / ( item['side']**2 )
        
        item['h_x'] = h_x_all[ h_all > 0 ]
        item['h']   = h_all  [ h_all > 0 ]

        item['thepath'] = re.match(".+\\\\([^\\\\]+\\\\[^\\\\]+)\\\\.+", str( d )).groups()[0]

        data.append( item )

    return np.array( data )

def remove_duplicates( data, key, discriminant ):
    choosen_idxs = {}
    for i, d in enumerate( data ):
        thekey = key( d )
        if( thekey in choosen_idxs ):
            if( discriminant( d ) > discriminant( data[i] ) ):
                choosen_idxs[ thekey ] = i
        else:
            choosen_idxs[ thekey ] = i

    return data[ list( choosen_idxs.values() ) ]

