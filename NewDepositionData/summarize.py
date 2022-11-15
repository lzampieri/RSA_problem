import sys
from glob import glob
import json
import os
from pathlib import Path
import numpy as np
import pandas as pd
import pyarrow

assert len(sys.argv) > 1
path = sys.argv[1]
to_remove = ['.','/','\\']
while( path[0] in to_remove ):
    path = path[1:]
while( path[-1] in to_remove ):
    path = path[:-1]

print("Summarizing in folder",path)

folders = list( glob( path + "/f_*/" ) )
print("Found", len(folders), "folders")

def load_file( filename ):
    text = filename.read_text()
    text = text.replace("nan","0")
    text = text.replace("False", "false" )
    text = text.replace("True", "true" )
    text = text.replace(",]", "]" )
    return json.loads( text )

def compute_in_decades( thearray, func, method = 'auto' ):
    num_sectors = 10
    tot_len = len( thearray )
    extremities = np.linspace( 0, tot_len, num_sectors, dtype = int )
    thearray = np.array( thearray )

    results = []
    for i in range( num_sectors - 1 ):
            
        data = thearray[ extremities[i] : extremities[i+1] ]
        
        results.append( func( data ) )
    
    result = func( thearray )

    if( method == 'ptp' ):
        std = np.ptp( results ) / 2
    elif( method == 'uniform' ):
        std = np.ptp( results ) / np.sqrt( 12 )
    elif( method == 'std' ):
        std = np.std( results, ddof = 1 )
    elif( method == 'auto' ):
        std = np.max( [ np.std( results, ddof = 1 ), np.ptp( results ) / np.sqrt( 12 ) ] )
    else:
        raise ValueError( f"{method} is not a recognized method. Use ptp, uniform, std or auto" )
    
    return result, std

data = []

for f in folders:

    d = Path( os.path.dirname( f ) )
    
    if( not (
        ( d / 'details.txt').exists() and
        ( d / 'deposition.txt').exists() ) ):
        print("Folder", f, "skipped")

    item = load_file( d / 'details.txt' )

    chunks = load_file( d / 'deposition.txt' )["occupation_history"]

    item['mean_v'],item['mean_s'] = compute_in_decades( chunks, lambda arr: np.mean( arr ) )
    item['std_v'],item['std_s'] = compute_in_decades( chunks, lambda arr: np.std( arr, ddof = 1 ) )

    data.append( item )

data = pd.DataFrame.from_records( data )

data.to_parquet( d / 'summary.parquet' )

print("Done!")