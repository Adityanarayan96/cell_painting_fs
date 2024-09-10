import logging

import pandas as pd
from pycytominer.operations import correlation_threshold, variance_threshold

from .metadata import find_feat_cols

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
import heapq
import json
import numpy as np

def select_features(dframe_path, feat_selected_path):
    '''Run feature selection'''
    dframe = pd.read_parquet(dframe_path)
    features = find_feat_cols(dframe.columns)
    low_variance = variance_threshold(dframe, features)
    features = [f for f in features if f not in low_variance]
    logger.info(f'{len(low_variance)} features removed by variance_threshold')
    high_corr = correlation_threshold(dframe, features, threshold=0.9)
    features = [f for f in features if f not in high_corr]
    logger.info(f'{len(high_corr)} features removed by correlation_threshold')

    dframe.drop(columns=low_variance + high_corr, inplace=True)

    dframe.reset_index(drop=True).to_parquet(feat_selected_path)

def top_intersecting_features(dframe_path, feat_selected_path_bio, feat_selected_path_batch, output_dframe_path , top_k_bio = 500, top_k_batch = 500):
    """
    Function that 1) selects top k features based on bio-metric and batch_alignment scores
    2) Selects features that are common between the two scores
    and 3) Saves the resulting dataframe to a parquet file

    Parameters
    ----------
    dframe_path : str
        Path to the input dataframe (usually the parquet file with imputed features)
    feat_selected_path : str
        Path to the dictionary with the selected features (_bio and _batch)
    output_dframe_path : str
        Path to the output dataframe (parquet file)
    top_k_bio : int
        number of features included from bio-metric score
    top_k_batch : int
        number of features included from batch_alignment score

    Future work wil update to select number of features you want to retain
    """
    dframe = pd.read_parquet(dframe_path) # Read the input dataframe
    # features = find_feat_cols(dframe.columns) # Get the feature columns
    with open(feat_selected_path_bio, 'r') as file:
        bio_dict = json.load(file)
    with open(feat_selected_path_batch, 'r') as file:
        batch_dict = json.load(file)

    bio_features = get_top_k_values(bio_dict, top_k_bio) # Get the top k bio-metric features
    batch_features = get_top_k_values(batch_dict, top_k_batch) # Get the top k batch_alignment features

    common_features =list(bio_features.keys() & batch_features.keys() & set(dframe.columns)) # Get the common features between the two scores
    meta_features = [f for f in dframe.columns if f.startswith('Metadata_')] # Get the metadata features
    
    common_features = common_features + meta_features # Add the metadata features to the common features

    logger.info(f'{len(common_features)} selected features are common between bio-metric and batch_alignment scores')
    dframe = dframe[common_features] # Select the common features
    # features_to_drop = [f for f in features if f not in common_features] # Get the features to drop
    # dframe.drop(columns=features_to_drop, inplace=True) # Drop the features that are not common between the two scores

    dframe.reset_index(drop=True).to_parquet(output_dframe_path)

def get_top_k_values(input_dict, k):
    # Get the top k keys with the largest values
    top_k_keys = heapq.nlargest(k, input_dict, key=input_dict.get)
    # Create a new dictionary with these top k keys
    top_k_dict = {key: input_dict[key] for key in top_k_keys}
    return top_k_dict

def add_random_features_to_top_k(top_k_features_path, k_1):
    # Load the existing top_k_features
    top_k_features = pd.read_parquet(top_k_features_path)

    # Generate k_1 random features
    random_features = pd.DataFrame(
        data=np.random.normal(0, 1, size=(len(top_k_features), k_1)),
        columns=[f'Random_Feature_{i}' for i in range(1, k_1 + 1)]
    )
    
    # Append the random features to the top_k_features
    top_k_features = pd.concat([top_k_features, random_features], axis=1)
    
    # Save the modified top_k_features back to top_k_features_path
    top_k_features.to_parquet(top_k_features_path, index=False)

    