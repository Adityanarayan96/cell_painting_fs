import sys
sys.path.append('/home/anravi2/cell_painting/Batch_effects_viewer/src')

from importlib import reload
import os
from metrics import scib as scib_metrics
# from methods import track_progress
# Add the custom package path
import warnings
warnings.filterwarnings("ignore")


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
    common_folder = "../../scenario_2"

    # Create the common folder if it doesn't exist
    if not os.path.exists(common_folder):
        os.makedirs(common_folder)

    # # File paths
    # parquet_path = "/home/anravi2/cell_painting/2023_Arevalo_BatchCorrection/outputs/scenario_1/raw.parquet" #Merged metadata
    # stats_path = os.path.join(common_folder, "test.parquet")
    # neg_stats_path = os.path.join(common_folder, "test_neg.parquet")
    # variant_feats_path = os.path.join(common_folder, "test_variant_feats.parquet")
    # normalize_path = os.path.join(common_folder, "test_mad.parquet")
    # outliers_path = os.path.join(common_folder, "test_outliers.parquet")
    # drop_path = os.path.join(common_folder, "test_drop.parquet")
    # clip_path = os.path.join(common_folder, "test_clip.parquet")
    # int_path = os.path.join(common_folder, "test_int.parquet")
    # feat_path = os.path.join(common_folder, "test_feat.parquet")
    # imputeknn_path = os.path.join(common_folder, "test_imputeknn.parquet")
    imputemedian_path = os.path.join(common_folder, "imputemedian.parquet")

    adata_path_feature = os.path.join(common_folder, "clusters_featurewise.h5ad")

    scib_metrics.cluster_featurewise(imputemedian_path, adata_path_feature)

    nmi_path = os.path.join(common_folder, "nmi_featurewise")
    label_key = "Metadata_JCP2022"
    # batch_key = "Metadata_Batch"

    scib_metrics.nmi_featurewise(adata_path_feature, label_key, nmi_path)
    

if __name__ == '__main__':
    main()