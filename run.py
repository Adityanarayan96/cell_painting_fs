import os

columns = [
    "Metadata_Source",
    "Metadata_Plate",
    "Metadata_Well",
    "Cells_AreaShape_Eccentricity",
    "Nuclei_AreaShape_Area",
]

columns_str = ','.join(columns)

experiments = [f"--metadata_source=source_6  --metadata_batch=L04  --metadata_plate=110000295627 --metadata_platetype=TARGET2 --columns={columns_str} ",
               ]

command_base = "python3 ./src/main_imageviewer.py {args} --git_clone_dir ../2023_Arevalo_BatchCorrection/inputs"

for i, experiment in enumerate(experiments):
    command = command_base.format(args=experiment)
    print(f"Executing command: {command}")
    os.system(command)