# Update the profile_formatter and loaddata_formatter using the command line arguments
import io
import os
import pandas as pd
import plotly.express as px
import plotly.io as pio
import pyarrow as pl
from config import *
from methods import *
from pickle import FALSE, TRUE
import requests
from io import BytesIO
from matplotlib import pyplot as plt
from matplotlib import image as mpimg
import boto3
from botocore import UNSIGNED
from botocore.config import Config

# plt.ion()
def main():
    #These are foramtters to obtain data from AWS servers
    profile_formatter = (
        "s3://cellpainting-gallery/cpg0016-jump/"
        "{Metadata_Source}/workspace/profiles/"
        "{Metadata_Batch}/{Metadata_Plate}/{Metadata_Plate}.parquet"
    )

    loaddata_formatter = (
        "s3://cellpainting-gallery/cpg0016-jump/"
        "{Metadata_Source}/workspace/load_data_csv/"
        "{Metadata_Batch}/{Metadata_Plate}/load_data_with_illum.parquet"
    )

    GIT_CLONE_DIR = args.git_clone_dir #Directory where metadata is stored
    file = pd.read_csv("/home/anravi2/cell_painting/Batch_effects_viewer/clicked_points.csv", dtype=str)
  
    ann_dframe = load_data_from_single_file(file, loaddata_formatter)
    print(ann_dframe.columns)
    show_images_single_file_test(ann_dframe, combine=True)
if __name__ == '__main__':
    main()
