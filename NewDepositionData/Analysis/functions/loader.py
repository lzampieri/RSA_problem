import numpy as np
import json
from pathlib import Path
import os
from glob import glob
import re
import pandas as pd
from . import stats

# Load single file
def load_file( filename ):
    text = filename.read_text()
    text = text.replace("nan","0")
    text = text.replace("False", "false" )
    text = text.replace("True", "true" )
    text = text.replace(",]", "]" )
    return json.loads( text )

def load_all_data():
    return load_data( "../*Analysis*/**/" )

def load_data( regex =  "../CloudVenetoAnalysisHuge/**/" ):
    data = []

    for file in glob( regex + "details.txt", recursive=True):
        item = {}
        d = Path( os.path.dirname( file ) )

        # Verify that everything exists
        # if( not (
        #     ( d / 'details.txt').exists() and
        #     ( d / 'deposition.txt').exists() and
        #     ( d / 'chunks.txt').exists() ) ):
        if( not (
            ( d / 'details.txt').exists() and
            ( d / 'deposition.txt').exists() ) ):
            continue

        # Verify that the item should not be excluded
        if( ( d / 'ignore.txt').exists() ):
            continue

        # Load data
        item.update( load_file( d / 'details.txt'     ) )
        item.update( load_file( d / 'deposition.txt'  ) )

        if( ( d / 'percolation.txt').exists() ):
            item.update( load_file( d / 'percolation.txt' ) )

        # with open( d / "chunks.txt" ) as file:
        #     lines = file.readlines()
        #     if( len( lines ) > 1 ):
        #         chunks = json.loads( "[" + lines[-1].replace(",]","]").replace("nan","0") + "]" )[-1]
        #     else:
        #         chunks = item['occupation_history']
        item['chunks'] = item['occupation_history'] # only for retrocompatibility

        # Correct parity
        if( item['side'] % 2 == 1 ):
            item['side'] = item['side'] - 1


        # Prepare histogram
        bins = np.arange( min( chunks ), max( chunks ) + 1 ) - 0.5
        h_all, _ = np.histogram( chunks, bins= bins, density=True )
        h_x_all = ( bins[:-1] + 0.5 ) / ( item['side']**2 )
        
        item['xs'] = np.array( chunks ) / ( item['side']**2 )
        item['h_x'] = h_x_all[ h_all > 0 ]
        item['h']   = h_all  [ h_all > 0 ]

        # print( str( d ) )

        match = re.match(".+\\\\([^\\\\]+\\\\[^\\\\]+)\\\\.+", str( d ))
        item['thepath'] = match.groups()[0] if match else str( d )

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

def compute_stds( data ):
    for d in data:
        d['jamming_frac_avg'] = stats.compute_in_decades(
            d['chunks'],
            lambda arr: np.mean( arr ) / d['side'] / d['side']
        )
        d['jamming_frac_std'] = stats.compute_in_decades(
            d['chunks'],
            lambda arr: np.std( arr, ddof = 1 ) / d['side'] / d['side']
        )

def filter( data, func ):
    data = np.array( data )
    return data[ np.vectorize( func )( data ) ]

def export( data, filename, columns = ['side','defects_frac','gamma','dep_polymers','runned_replicas','chunks'], renames = {} ):
    df = pd.DataFrame.from_records( data )
    df = df[columns]
    renames_int = { 'runned_replicas': 'total_replicas', 'chunks': 'jamming_occupations' }
    renames_int.update( renames )
    df.rename( columns=renames_int, inplace=True )
    df.to_json( filename, orient='records', indent = 2 )
    print( len( df ), " rows exported to ", filename )