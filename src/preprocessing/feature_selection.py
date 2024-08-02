# Description: This file contains functions for feature selection.
import plotly.express as px
import plotly.io as pio
import pyarrow as pl
from matplotlib import pyplot as plt
from matplotlib import image as mpimg

def feature_metrics_list(dframe_path):
    """
    This function creates a list of metrics calculated for single feature combinations 
    """

    



def feature_metric_calculation(dframe_path, indice, metric):
    """
    This function performs feature selection on the data.
    data: A pandas dataframe.
    indice: indice of feature to calculate metrics on
    metric to calculate

    """

    # Calculate the metric
    metric_value = metrics.metric(dframe_path)

    return metric_value