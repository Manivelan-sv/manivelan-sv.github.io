'''
Dynamical instabilities cause extreme events in a theoretical Brusselator model.
Manivelan, S. V., Sabarathinam, S., Thamilmaran, K., & Manimehan, I.
Chaos, Solitons & Fractals, 180, 114582. (https://doi.org/10.1016/j.chaos.2024.114582)
'''

import numpy as np

data = []  # Initialize an empty list to store data

def create_dat_file(input_file, output_file):
    # Open the input file for reading
    with open(input_file, 'r') as f:
        # Read each line from the input file
        for line in f:
            # Split the line into two values and convert them to floats
            t, x = map(float, line.split())
            # Append the second value (x) to the data list
            data.append(x)

    # Create pairs of consecutive values from the data list
    pairs = zip(data, data[1:])

    # Open the output file for writing
    with open(output_file, 'w') as f:
        # Write each pair of values into the output file in tab-separated format
        for pair in pairs:
            f.write(f"{pair[0]}\t{pair[1]}\n")
    
    # Print a message indicating that the .dat file has been created successfully
    print(f".dat file '{output_file}' created successfully.")

# Usage example
input_file = '1.1125.dat'      # Replace with the actual input file name
output_file = '1.1125_rm.dat'  # Replace with the desired output file name

# Call the function to create the .dat file
create_dat_file(input_file, output_file)
