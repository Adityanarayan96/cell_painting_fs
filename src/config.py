import argparse

# Define the command line arguments
parser = argparse.ArgumentParser(description='Batch Effects Viewer Configuration')
# Directory with data
parser.add_argument('--git_clone_dir', type=str, help='Git clone directory') 
# Metadata source, leave empty if no specific source
parser.add_argument('--metadata_source', type=str, help='Metadata_Source') 
# Metadata batch, leave empty if no specific batch
parser.add_argument('--metadata_batch', type=str, help='Metadata_Batch')
# Metadata plate, leave empty if no specific plate
parser.add_argument('--metadata_plate', type=str, help='Metadata_Plate')
# Metadata well, leave empty if no specific well
parser.add_argument('--metadata_platetype', type=str, help='Metadata_PlateType')
# Columns to load
parser.add_argument('--columns', type=lambda s: [str(item) for item in s.split(',')], help='Columns to load or None')
args = parser.parse_args() #Parses the command line arguments

if args.columns == "None":
    args.columns = None
