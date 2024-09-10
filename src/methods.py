import os
import sys

current_dir = os.path.dirname(os.path.abspath(__file__))
if current_dir not in sys.path:
    sys.path.append(current_dir)

# import io
import random
import os
import pandas as pd
import plotly.express as px
import plotly.io as pio
import pyarrow as pl
# import plotly.express as px
import plotly.graph_objects as go
# from config import args
from pickle import FALSE, TRUE
# import requests
# from io import BytesIO
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
import scanpy as sc
import anndata as ad
import plotly.graph_objects as go
from ipywidgets import widgets
from IPython.display import display


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




# Method to display an interactive UMAP plot with metadata labels
def display_interactive_umap(adata_path, clicked_points_path):
    adata = ad.read_h5ad(adata_path)
    meta = adata.obs

    # Compute UMAP embedding
    sc.tl.umap(adata)
    umap_coords = adata.obsm['X_umap']

    # Initialize the click data list and DataFrame
    click_data = []
    click_df = pd.DataFrame(columns=['Metadata_Source', 'Metadata_Batch', 'Metadata_Plate', 'Metadata_Well'])

    # Dropdown widget for selecting the metadata column
    metadata_cols = sorted(meta.columns.tolist())  # Sort columns lexicographically
    dropdown = widgets.Dropdown(
        options=metadata_cols,
        value=metadata_cols[0],
        description='Label:',
    )

    # Function to create the scatter plot with a specific metadata column for labeling and color coding
    def create_scatter_plot(label_col):
        # Map the metadata column to colors using Plotly's built-in function
        color_sequence = px.colors.qualitative.Plotly
        unique_values = sorted(meta[label_col].unique())  # Sort values lexicographically
        color_map = {val: color_sequence[i % len(color_sequence)] for i, val in enumerate(unique_values)}

        fig = go.FigureWidget()

        for value in unique_values:
            indices = meta[label_col] == value
            scatter = go.Scatter(
                x=umap_coords[indices, 0] * 5, 
                y=umap_coords[indices, 1],
                mode='markers',
                marker=dict(color=color_map[value]),
                name=value,
                text=meta[indices][label_col],
                customdata=meta[indices][['Metadata_Source', 'Metadata_Batch', 'Metadata_Plate', 'Metadata_Well']],
                visible=True
            )
            fig.add_trace(scatter)

        fig.update_layout(
            title="UMAP Scatter Plot",
            updatemenus=[
                {
                    "buttons": [
                        {
                            "label": value,
                            "method": "update",
                            "args": [
                                {"visible": [val == value for val in unique_values]},
                                {"title": f"UMAP Scatter Plot - {value}"}
                            ]
                        }
                        for value in unique_values
                    ],
                    "direction": "down",
                    "showactive": True
                }
            ]
        )

        # Add click event handling for all traces
        for trace in fig.data:
            trace.on_click(handle_click)

        return fig

    # Click event handling function
    def handle_click(trace, points, state):
        nonlocal click_data, click_df
        if points.point_inds:
            ind = points.point_inds[0]
            click_info = {
                'Metadata_Source': trace.customdata[ind][0],
                'Metadata_Batch': trace.customdata[ind][1],
                'Metadata_Plate': trace.customdata[ind][2],
                'Metadata_Well': trace.customdata[ind][3]
            }
            click_data.append(click_info)

            # Convert the click_info to a DataFrame and concatenate
            click_info_df = pd.DataFrame([click_info])
            click_df = pd.concat([click_df, click_info_df], ignore_index=True)
            click_df.to_csv(clicked_points_path, index=False)

            print("Click registered:", click_info)

    # Function to update the plot based on dropdown selection
    def update_plot(change):
        fig = create_scatter_plot(change.new)
        display(fig)

    dropdown.observe(update_plot, names='value')

    # Display the dropdown and the initial plot
    display(dropdown)
    initial_fig = create_scatter_plot(dropdown.value)
    display(initial_fig)


def create_manhattan_plot(grouping, save_path= "manhattan.html", title = 'Manhattan Plot with Scatter'): #Input is a grouped dictionary
    # Map each unique group name to a numeric value
    unique_groups = list(grouping.keys())
    group_to_numeric = {group: i for i, group in enumerate(unique_groups)}

    # Prepare data for plotting with slight scatter in x-axis values
    x_values = []
    y_values = []
    hover_texts = []
    marker_colors = []

    colors = ['red', 'blue', 'green', 'orange', 'purple', 'brown', 'pink', 'gray', 'cyan', 'magenta']
    group_to_color = {group: colors[i % len(colors)] for i, group in enumerate(unique_groups)}

    # Loop through each group and collect the values, adding scatter to x-values
    for group, items in grouping.items():
        for key, value in items:
            # Convert the group to a numeric value and add random scatter
            numeric_x = group_to_numeric[group] + random.uniform(-0.4, 0.4)
            x_values.append(numeric_x)
            y_values.append(value)
            hover_texts.append(key)
            marker_colors.append(group_to_color[group])  

    # Create the Manhattan plot
    fig = go.Figure()

    fig.add_trace(go.Scatter(
        x=x_values,
        y=y_values,
        mode='markers',
        marker=dict(size=10, color=marker_colors),
        text=hover_texts,
        hoverinfo='text',
    ))

    # Set plot title and labels
    fig.update_layout(
        title=title,
        xaxis=dict(
            tickmode='array',
            tickvals=list(group_to_numeric.values()),
            ticktext=list(group_to_numeric.keys())
        ),
        yaxis_title='Values',
        showlegend=False
    )
    #Save the figure
    fig.write_html(save_path)
    # Display the plot
    fig.show()