import os
import sys

current_dir = os.path.dirname(os.path.abspath(__file__))
if current_dir not in sys.path:
    sys.path.append(current_dir)

import io
import os
import pandas as pd
import plotly.express as px
import plotly.io as pio
import pyarrow as pl
# import plotly.express as px
import plotly.graph_objects as go
from config import args
from pickle import FALSE, TRUE
import requests
from io import BytesIO
from matplotlib import pyplot as plt
from matplotlib import image as mpimg
import boto3
from botocore import UNSIGNED
from botocore.config import Config
import anndata as ad
import scanpy as sc
from PIL import Image
import numpy as np
from matplotlib.colors import Normalize
from skimage.color import rgb2gray
import skimage.io
import os
import boto3
from botocore import UNSIGNED
from botocore.client import Config
import matplotlib.image as mpimg
import numpy as np
from io import BytesIO
import skimage.io

def show_images_single_file_test(file, combine=False):
    channels = ['OrigDNA', 'OrigAGP', 'OrigER', 'OrigMito', 'OrigRNA']
    channel_colors = {
        'OrigDNA': (0, 0, 1),  # Blue
        'OrigAGP': (0, 0.3, 0),  # Dark Green
        'OrigER': (0, 0.4, 0),  # Medium Green
        'OrigRNA': (0, 0.5, 0),  # Light Green
        'OrigMito': (1, 0, 0)  # Red
    }
    s3_client = boto3.client("s3", config=Config(signature_version=UNSIGNED))

    for _, row in file.iterrows():
        directory_path = f"images/{row['Metadata_Source']}/{row['Metadata_Plate']}/{row['Metadata_Well']}/{row['Metadata_Site']}"
        os.makedirs(directory_path, exist_ok=True)
        combined_image = None

        for channel in channels:
            path_name_column = f'PathName_{channel}'
            file_name_column = f'FileName_{channel}'
            image_url = os.path.join(row[path_name_column], row[file_name_column])

            bucket = image_url.split("/")[2]
            key = "/".join(image_url.split("/")[3:])
            
            response = s3_client.get_object(Bucket=bucket, Key=key)
            image = mpimg.imread(BytesIO(response["Body"].read()), format="tiff") 
            image = image / image.max()
            
            # Save individual channel images
            if combine is False:
                output_path = os.path.join(directory_path, f"{channel}_{row['Metadata_Well']}_{row['Metadata_Site']}.tiff")
                mpimg.imsave(output_path, image, format="tiff")

            if combine is True:
                if combined_image is None:
                    combined_image = np.zeros((*image.shape, 3))

                color = channel_colors[channel]
                for i in range(3):
                    combined_image[:, :, i] += image * color[i]
        
        if combine is True:
            combined_image = combined_image/ combined_image.max() # Normalize to avoid overflow
            combined_output_path = os.path.join(directory_path, f"combined_{row['Metadata_Well']}_{row['Metadata_Site']}.tiff")
            mpimg.imsave(combined_output_path, combined_image)

def show_images_single_file(file, combine = False):
    channels = ['OrigDNA', 'OrigAGP', 'OrigER', 'OrigMito', 'OrigRNA']
    s3_client = boto3.client("s3", config=Config(signature_version=UNSIGNED))

    for _, row in file.iterrows():
        directory_path = f"images/{row['Metadata_Source']}/{row['Metadata_Plate']}/{row['Metadata_Well']}/{row['Metadata_Site']}"
        os.makedirs(directory_path, exist_ok=True)
        for channel in channels:
            path_name_column = f'PathName_{channel}'
            file_name_column = f'FileName_{channel}'
            image_url = os.path.join(row[path_name_column], row[file_name_column])

            bucket = image_url.split("/")[2]
            key = "/".join(image_url.split("/")[3:])
            
            response = s3_client.get_object(Bucket=bucket, Key=key)
            image = mpimg.imread(BytesIO(response["Body"].read()), format="tiff") 
            image         
            combined_image = []
            # Save individual channel images
            if combine is False:
                output_path = os.path.join(directory_path, f"{channel}_{row['Metadata_Well']}_{row['Metadata_Site']}.tiff")
                mpimg.imsave(output_path, image, format="tiff")

            if combine is True:
                combined_image.append(image)
            
        if combine is True:
            combined_image = np.array(combined_image)
            # combined_image = combined_image/combined_image.max()
            combined_output_path = os.path.join(directory_path, f"combined_{row['Metadata_Well']}_{row['Metadata_Site']}.tiff")
            skimage.io.imsave(combined_output_path, combined_image)

# Method to extract data based on the command line arguments
def get_by_plates(plates, num_of_samples=None):
    if plates.empty:  # Check if the 'plates' variable is empty
        raise ValueError("The plates variable is empty")
    if num_of_samples is not None:
        return plates.sample(num_of_samples)
    else:
        return plates

#Method to get data from a single file
def load_data_from_single_file(file, formatter=None, columns=None):
    """
    Use appropriate profiller to get either well level profiles or load data
    """
    if formatter is None:
        raise ValueError("formatter is required")
    dframes = []
    # filtered_data = pd.merge(merged_data, file, on=["Metadata_Source", "Metadata_Plate", "Metadata_Well"], how='inner')
 # Iterate over each row in filtered_data
    for _, row in file.iterrows():
        # Generate the S3 path from the row data
        s3_path = formatter.format(**row.to_dict())
        
        # Read the parquet file from S3
        df_parquet = pd.read_parquet(s3_path, storage_options={"anon": True}, columns=columns)
        
        # Perform an inner join between df_parquet and filtered_data using the keys
        # This assumes 'Metadata_Source', 'Metadata_Plate', and 'Metadata_Well' are present in both DataFrames
        # and can uniquely identify rows for the purpose of filtering
        df_filtered = pd.merge(df_parquet, file[['Metadata_Source', 'Metadata_Plate', 'Metadata_Well']],
                            on=['Metadata_Source', 'Metadata_Plate', 'Metadata_Well'], how='inner')
        
        # Append the filtered DataFrame to the list
        dframes.append(df_filtered)
    return pd.concat(dframes)
#Method to load well level profiles with features, also annotates metadata with well level profiles
def get_well_level_profiles(source ,filtered_plate, compound, wells , profile_formatter=None, columns=None, num_samples = None):
    if source == '':
        raise ValueError("The source variable is empty")
    if filtered_plate.empty:
        raise ValueError("The filtered_df variable is empty")
    if profile_formatter is None:
        raise ValueError("profile_formatter is required")
    if compound.empty:
        raise ValueError("The compound variable is empty")
    if wells.empty:
        raise ValueError("The wells variable is empty")
    dframes = load_data(source, filtered_plate, formatter=profile_formatter, columns=columns, num_samples=num_samples)
    # Join features from metadata
    metadata = compound.merge(wells, on="Metadata_JCP2022") #JCP2022 is the unique identifier
    ann_dframe = metadata.merge(
    dframes, on=["Metadata_Source", "Metadata_Plate", "Metadata_Well"]
    )
    return ann_dframe

#load data gets metadata to load images
def load_data(source ,filtered_plate, formatter=None, columns=None, num_samples = None):
    if source == '':
        raise ValueError("The source variable is empty")
    if filtered_plate.empty:
        raise ValueError("The filtered_df variable is empty")
    if formatter is None:
        raise ValueError("formatter is required")
    dframes =[]
    filtered_plate = filtered_plate[filtered_plate['Metadata_Source'] == source]
    filtered_plate = filtered_plate.sample(num_samples) if num_samples is not None else filtered_plate
    # Reads and stores parquet files (from aws server) based on filtered_plate (plate that you want wells from)
    for _, row in filtered_plate.iterrows():
        s3_path = formatter.format(**row.to_dict())
        dframes.append(
            pd.read_parquet(s3_path, storage_options={"anon": True}, columns=columns) #Reads parquet files anonymously
        )
    return pd.concat(dframes)

#Plot 2-D features indexed by metadata
def plot_features(data=None, x_feature="Cells_AreaShape_Eccentricity", 
                  y_feature="Nuclei_AreaShape_Area",
                  color="Metadata_Source", hover_name="Metadata_JCP2022", 
                  hover_data=["Metadata_InChIKey"],):
    if data is None:
        raise ValueError("The data variable is empty")
    fig = px.scatter(data, x=x_feature, y=y_feature, color=color, 
                     hover_name=hover_name, hover_data=hover_data)
    fig.show()
    # pio.write_html(fig, file="scatter_plot.html", auto_open=True)

#Method to display all fields of images and combined channel images from single metadata
# def show_images_single_file(file, combined=True):
#     channels = ['OrigDNA', 'OrigAGP', 'OrigER', 'OrigMito', 'OrigRNA']
#     color_map = {
#         'OrigDNA': 'Reds',   # Red colormap
#         'OrigAGP': 'Greens', # Green colormap
#         'OrigER': 'Blues',   # Blue colormap
#         'OrigMito': 'YlOrBr',# Yellow-orange-brown colormap
#         'OrigRNA': 'PuBuGn'  # Purple-blue-green colormap
#     }

#     s3_client = boto3.client("s3", config=Config(signature_version=UNSIGNED))

#     for _, row in file.iterrows():
#         directory_path = f"images/{row['Metadata_Source']}/{row['Metadata_Plate']}/{row['Metadata_Well']}/{row['Metadata_Site']}"
#         os.makedirs(directory_path, exist_ok=True)
        
#         if combined:
#             combined_image = None

#         for channel in channels:
#             path_name_column = f'PathName_{channel}'
#             file_name_column = f'FileName_{channel}'
#             image_url = os.path.join(row[path_name_column], row[file_name_column])

#             bucket = image_url.split("/")[2]
#             key = "/".join(image_url.split("/")[3:])
            
#             response = s3_client.get_object(Bucket=bucket, Key=key)
#             image = Image.open(BytesIO(response["Body"].read())).convert("L")
#             image = np.array(image)
#             if combined:
#                 # Normalize the image and add it to the combined image   
#                 image_normalized = image / np.max(image)
#                 if combined_image is None:
#                     combined_image = [image_normalized]
#                 else:
#                     combined_image = np.concatenate((combined_image, [image_normalized]), axis=2)
#             else:
#                 # Save individual channel images
#                 output_path = os.path.join(directory_path, f"{channel}_{row['Metadata_Well']}_{row['Metadata_Site']}.tiff")
#                 Image.fromarray(image).save(output_path, format="TIFF")

#         if combined and combined_image is not None:
#             combined_output_path = os.path.join(directory_path, f"combined_{row['Metadata_Well']}_{row['Metadata_Site']}.tiff")
#             Image.fromarray(combined_image).save(combined_output_path, format="TIFF")



#Method to display images corresponding to meta_data.
def show_images(linked_df, channel = None, iloc = 0): #Need to modify and add functionality
    """
    
    """
    # Assuming linked_df is your DataFrame and iloc is the index you're interested in
    channels = ['OrigDNA', 'OrigAGP', 'OrigER', 'OrigMito', 'OrigRNA']
    if channel is not None:
        path_name_column = f'PathName_{channel}'
        file_name_column = f'FileName_{channel}'
        image_url = os.path.join(
        linked_df.iloc[iloc][path_name_column], linked_df.iloc[iloc][path_name_column]
    )
        s3_client = boto3.client("s3", config=Config(signature_version=UNSIGNED)) #Configuring boto3 client
        response = s3_client.get_object(
        Bucket=image_url.split("/")[2], Key="/".join(image_url.split("/")[3:])
        )
        image = mpimg.imread(BytesIO(response["Body"].read()), format="tiff")

        # fig = plt.imshow(image)
        # print(image_url)
        # fig.show()

    else:
    # Create a figure with 1 row and 5 columns
        fig, axs = plt.subplots(1, 5, figsize=(20, 4))
        for i, channel in enumerate(channels):
            path_name_column = f'PathName_{channel}'
            file_name_column = f'FileName_{channel}'
            image_url = os.path.join(
                linked_df.iloc[iloc][path_name_column], linked_df.iloc[iloc][file_name_column]
            )

            s3_client = boto3.client("s3", config=Config(signature_version=UNSIGNED))
            response = s3_client.get_object(
                Bucket=image_url.split("/")[2], Key="/".join(image_url.split("/")[3:])
            )

            image = mpimg.imread(BytesIO(response["Body"].read()), format='tiff')

            axs[i].imshow(image)
            axs[i].set_title(channel)
            axs[i].axis('off')  # Hide axis

    # Save the plot to a file
    plt.savefig('output_figure.tiff', dpi=300, format='tiff')

    # # Display the saved plot
    # Image('output_figure.png')
    
    