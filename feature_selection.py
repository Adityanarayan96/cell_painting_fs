import io
import os
import pandas as pd
import plotly.express as px
import plotly.io as pio
import pyarrow as pl
from src.config import *
from src.methods import *
from pickle import FALSE, TRUE
import requests
from io import BytesIO
from matplotlib import pyplot as plt
from matplotlib import image as mpimg
import boto3
from botocore import UNSIGNED
from botocore.config import Config


def main():
    #These are formatters to obtain data from AWS servers
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

if __name__ == '__main__':
    main()