import os
import csv

# Define the path to the output directory
output_dir = 'output'
output_csv = 'extracted_data.csv'

# Initialize a list to store the extracted data
data_list = []

# Traverse the directory
for subdir, _, files in os.walk(output_dir):
    for file in files:
        if file == 'CL_CD.txt':
            file_path = os.path.join(subdir, file)
            with open(file_path, 'r') as f:
                # Read lines from the file and extract required data
                lines = f.readlines()
                mach = float(lines[0].split(':')[1].strip())
                aoa = float(lines[1].split(':')[1].strip())
                cl = float(lines[2].split(':')[1].strip())
                cd = float(lines[3].split(':')[1].strip())
                
                # Append the extracted data to the list
                data_list.append([aoa, mach, cd, cl])

# Write the data to a CSV file
with open(output_csv, 'w', newline='') as csvfile:
    csvwriter = csv.writer(csvfile)
    # Write the header
    csvwriter.writerow(['AOA', 'MA', 'CD', 'CL'])
    # Write the data
    csvwriter.writerows(data_list)

print(f'Data has been successfully extracted and saved to {output_csv}')

