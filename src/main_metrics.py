import sys
sys.path.append('/home/anravi2/cell_painting/Batch_effects_viewer/src')

from importlib import reload
import os
from metrics import scib as scib_metrics
# from methods import track_progress
# Add the custom package path
import warnings
warnings.filterwarnings("ignore")

import anndata as ad
import pandas as pd

# from src.metrics import scib as scib_metrics


# Import and reload custom modules
import preprocessing.io as io
import preprocessing.stats as stats
import preprocessing.normalize as normalize
import preprocessing.outliers as outliers
import preprocessing.transform as transform
import preprocessing.feature_selection as feature_selection


reload(io)
reload(stats)
reload(normalize)
reload(outliers)
reload(transform)
reload(feature_selection)

# Import specific functions if needed
from preprocessing.io import *       # Assuming io.py contains necessary functions
from preprocessing.stats import *    # Assuming stats.py contains necessary functions
from preprocessing.normalize import * # Assuming normalize.py contains necessary functions
from preprocessing.outliers import *  # Assuming outliers.py contains necessary functions
from preprocessing.transform import * # Assuming transform.py contains necessary functions
from preprocessing.feature_selection import * # Assuming feature_selection.py contains necessary functions

def main():
    # File paths
    # Define the common folder path
    common_folder = "../../scenario_4"

    # Create the common folder if it doesn't exist
    if not os.path.exists(common_folder):
        os.makedirs(common_folder)

    # File paths
    parquet_path = "/home/anravi2/cell_painting/2023_Arevalo_BatchCorrection/outputs/scenario_4/raw.parquet" #Merged metadata
    stats_path = os.path.join(common_folder, "test.parquet")
    neg_stats_path = os.path.join(common_folder, "test_neg.parquet")
    variant_feats_path = os.path.join(common_folder, "test_variant_feats.parquet")
    normalize_path = os.path.join(common_folder, "test_mad.parquet")
    outliers_path = os.path.join(common_folder, "test_outliers.parquet")
    drop_path = os.path.join(common_folder, "test_drop.parquet")
    clip_path = os.path.join(common_folder, "test_clip.parquet")
    int_path = os.path.join(common_folder, "test_int.parquet")
    feat_path = os.path.join(common_folder, "test_feat.parquet")
    imputeknn_path = os.path.join(common_folder, "test_imputeknn.parquet")
    imputemedian_path = os.path.join(common_folder, "imputemedian.parquet")

    # adata_path_feature = os.path.join(common_folder, "clusters_featurewise.h5ad")

    # Function calls
    compute_negcon_stats(parquet_path, neg_stats_path)
    select_variant_features(parquet_path, neg_stats_path, variant_feats_path, union=False)
    mad(variant_feats_path, neg_stats_path, normalize_path)
    compute_stats(normalize_path, stats_path)
    iqr(100.0, normalize_path, stats_path, outliers_path)

    # Additional function calls based on the rules
    drop_cols(normalize_path, outliers_path, drop_path)
    clip_cols(drop_path, outliers_path, 500.0, clip_path)
    rank_int(clip_path, int_path)
    # select_features(int_path, feat_path)
    impute_knn(int_path, outliers_path, imputeknn_path)
    impute_median(imputeknn_path, outliers_path, imputemedian_path)
    # scib_metrics.cluster_featurewise(imputemedian_path, adata_path_feature)

    # nmi_path = os.path.join(common_folder, "nmi_featurewise")
    label_key = "Metadata_JCP2022"
    batch_key = "Metadata_Batch"
    silhoutte_label_single_path = os.path.join(common_folder, "test_silhoutte_single_label.bin")
    silhoutte_batch_single_path = os.path.join(common_folder, "test_silhoutte_single_batch.bin")

    scib_metrics.asw_single_features(imputemedian_path, label_key, silhoutte_label_single_path)
    scib_metrics.silhouette_batch_single_features(imputemedian_path, label_key, batch_key, silhoutte_batch_single_path)
    # scib_metrics.nmi_featurewise(adata_path_feature, label_key, nmi_path)
    

if __name__ == '__main__':
    main()