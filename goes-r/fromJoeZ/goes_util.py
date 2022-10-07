'''
goes_util.py
Utilities for GOES processing
'''

import pandas as pd

# create goes constants dataframe
def create_goes_constants_dataframe():
    arrays = [['goes09'],
         ['4']]
    tuples = list(zip(*arrays))
    columns = ['nu','c1','c2','scaling_gain','scaling_bias','a','b']
    index = pd.MultiIndex.from_tuples(tuples, names=['satellite', 'channel'])
    goes_df = pd.DataFrame(index=index,columns=columns)

    goes09_4 = [934.55, 1.191066e-5, 1.438833, 5.2285, 15.6854, -0.377608, 1.001284]
    goes_df.loc['goes09','4'] = goes09_4

    return goes_df

def get_goes_constants(goesNumber,channel):
    channel = str(channel)
    goes_df = create_goes_constants_dataframe()
    
    goes_constants = {}
    goes_constants['channel'] = channel
    goes_constants['nu'] = goes_df.loc[goesNumber,channel]['nu']
    goes_constants['c1'] = goes_df.loc[goesNumber,channel]['c1']
    goes_constants['c2'] = goes_df.loc[goesNumber,channel]['c2']
    goes_constants['scaling_gain'] = goes_df.loc[goesNumber,channel]['scaling_gain']
    goes_constants['scaling_bias'] = goes_df.loc[goesNumber,channel]['scaling_bias']
    goes_constants['a'] = goes_df.loc[goesNumber,channel]['a']
    goes_constants['b'] = goes_df.loc[goesNumber,channel]['b']
    return goes_constants
    

