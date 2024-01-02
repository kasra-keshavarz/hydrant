"""
Common functions for manupulating the gis outfput from other sources
"""

import pandas as pd
import numpy as np
from typing import Dict, Union
from itertools import product
import sys

def manipulating_fractions (df,
                            df_mapping: Dict[str, str] = {'prefix': 'Frac_'}, # {'id':'ID', 'prefix': 'Frac_'}
                            action = 'majority',
                            pass_maximum_fract_to_zero_rows = True,
                            minimum_value = 0):
    
    # Set the ID column to index for the DataFrame
    ID = df_mapping.get('id')
    if ID is not None:
        df.index = df[ID].astype(int)
        df.drop(columns=[ID], inplace=True)
        df = df.sort_index()
    
    # Extract columns that start with a prefix
    prefix = df_mapping.get('prefix')
    if prefix is not None:
        prefix_cols = [col for col in df.columns if col.startswith(prefix)]
    else:
        prefix_cols = df.columns
        
    # check the negative values for the prefix_cols
    all_positive = (df[prefix_cols] >= 0).all().all()
    if not all_positive:
        sys.exit("There are negative values in the Fractions")
    
    # find global maximume column
    df_sum = df.sum(axis=0)
    all_max_column = df_sum[prefix_cols].idxmax()
    
    # iterate over the row for each case of majority or normalize
    if action.lower() == 'majority':
        for index, row in df.iterrows():
            if sum (row[prefix_cols]) > 0:
                max_column = row[prefix_cols].idxmax()
                # print(max_column)
                df.loc[index,prefix_cols] = 0
                df.loc[index,max_column] = 1
            elif sum (row[prefix_cols]) == 0:
                if pass_maximum_fract_to_zero_rows:
                    df.loc[index,prefix_cols] = 0
                    df.loc[index,all_max_column] = 1
            else:
                sys.exit('sum of the rows are negative; fraction cannot be negative')
    elif action.lower() == 'normalize':
        for index, row in df.iterrows():
            if sum (row[prefix_cols]) > 0:
                if minimum_value is not None:
                    # check minimume to the maximume values
                    if minimum_value > df.loc[index,prefix_cols].max():
                        sys.exit("the maximume fraction of land cover is smaller than \
                        minimume_values provided for row/index "+str(index))
                    else:
                        df.loc[index, prefix_cols] = df.loc[index, prefix_cols].mask(df.loc[index, prefix_cols] < minimum_value, 0)
                df.loc[index,prefix_cols] = df.loc[index,prefix_cols]/sum (df.loc[index, prefix_cols])
            elif sum (row[prefix_cols]) == 0:
                if pass_maximum_fract_to_zero_rows:
                    df.loc[index,prefix_cols] = 0
                    df.loc[index,all_max_column] = 1
            else:
                sys.exit('sum of the rows are negative; fraction cannot be negative')
    else:
        sys.exit('action is not majority or normalize')
        
    df_modified = df.copy()
    
    # return
    return df_modified

def intersect_df(*dfs: pd.DataFrame,
                 df_mappings: Union[Dict[str, Dict[str, str]], None] = None,
                 remove_zero_combinations: bool = True):
    
    # turn input into list
    dfs = list(dfs)
    
    # Check if mappings are provided
    if df_mappings is None:
        df_mappings = {}

    # Loop through each DataFrame
    for idx, df in enumerate(dfs, start=1):
        mapping = df_mappings.get(f'df{idx}', {})
        ID = mapping.get('id')
        prefix = mapping.get('prefix', '')
        data_name = mapping.get('data_name', f'Data_{idx}')
        
        # Set the ID column to index for the DataFrame
        if ID is not None:
            df.index = df[ID].astype(int)
            df.drop(columns=[ID], inplace=True)
            df = df.sort_index()
        
        # Select columns based on prefix and replace them with data_name
        prefix_cols = [col for col in df.columns if col.startswith(prefix)]
        df = df[prefix_cols] # keep the columns with prefix
        for col in df.columns:
            new_col_name = col.replace(prefix, data_name)
            df[new_col_name] = df[col]
            df = df.drop(columns=[col])
            
        # update the dataframe
        dfs[idx - 1] = df
        
    # Compare the indexes
    index_comparison = all(dfs[i].index.equals(dfs[i + 1].index) for i in range(len(dfs) - 1))
    if index_comparison:
        print("The indexes of all DataFrames are exactly the same with the same order.")
    else:
        sys.exit("The indexes of DataFrames are different or in a different order.")
    
    # Generate combinations and calculate products
    combinations = list(product(*(df.columns for df in dfs)))
    result = pd.DataFrame(columns=[f'{" ".join(cols)}' for cols in combinations], \
                          index=dfs[0].index)
    combination_list = []
    
    # loop over combination
    for cols in combinations:
        col_product = dfs[0][cols[0]] #initialize the product
        for i in range(1, len(dfs)):
            col_product = col_product * dfs[i][cols[i]] # multiply the product
        # update the results
        result[f'{" ".join(cols)}'] = col_product
        
        # keep the positive or existing combinations
        if (sum(col_product)>0):
            combination_list.append(' '.join(cols))
        elif (sum(col_product)==0):
            if not remove_zero_combinations:
                combination_list.append(' '.join(cols))
        else:
            sys.exit("sum of columns cannot be negative")
    
    # total number of non zero combinations:
    print('total number of non zero combinations: ',len(combination_list)) 
    
    # Report combinations with non-zero values
    report = pd.DataFrame(combination_list, columns=['Combinations'])
    # Split combinations into separate columns based on DataFrame columns
    for idx, df in enumerate(dfs, start=0):
        mapping = df_mappings.get(f'df{idx+1}', {})
        data_name = mapping.get('data_name', f'Data_{idx+1}')
        report[data_name] = report['Combinations'].apply(lambda x: x.split()[idx])
        report[data_name] = report[data_name].str.extract(r'(\d+)')  # f'{chr(65 + i)}'
    report ['comb'] = np.arange(len(report))+1
    
    # remove the zero combination from result dataframe
    if remove_zero_combinations:
        zeros_columns = result.columns[(result == 0).all()]  # Get columns with all zeros
        # Drop columns that consist entirely of zeros
        result = result.drop(zeros_columns, axis=1)
        
    # rename the result header to 1 to n
    for i, col in enumerate(result.columns, start=1):
        new_col_name = f'comb_{i:04d}'  # Using f-string to format column names
        result.rename(columns={col: new_col_name}, inplace=True)
    
    return result , report