"""
Author: Kerry McGowan, 2022
"""

# Import required packages
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
# Magic command, display matplotlib plots inline
%matplotlib inline
import seaborn as sns
import numpy as np
import plotly as py
from plotly import graph_objects as go
from collections import Counter

# Directory location of FASTQ files
fastq_directory_R1 = './S1_S1_R1_001.fastq'
fastq_directory_R2 = './S1_S1_R2_001.fastq'

#### PART 1: Read in FASTQ files, extract sequencing reads into a data frame
# Define a function
def read_fastq_F(R1_fastq, R2_fastq, nrows):
    # Read in both FASTQ files
    R1_file = pd.read_csv(R1_fastq, header = None)
    R2_file = pd.read_csv(R2_fastq, header = None)
    
    # Extract sequencing reads (every 4th row)
    R1_sequences = R1_file.iloc[1::4,:]
    R2_sequences = R2_file.iloc[1::4,:]
    R1_sequences.columns = ['read1']
    R2_sequences.columns = ['read2']
    
    # Output a data frame with two columns containing each paired sequencing reads
    df = pd.concat([R1_sequences, R2_sequences], axis = 1)
    df.reset_index(drop = True, inplace = True)
    return(df.head(nrows))

# Apply function to the first 100,000 sequencing reads and save output as df1
df1 = read_fastq_F(fastq_directory_R1, fastq_directory_R2, 100000)
df1.to_csv('df1.csv')

#### PART 2: Extract barcodes at the 5' end of reads into a data frame
# Extract the 15 character barcodes from R1 and R2
read1 = pd.DataFrame(df1['read1'].str.slice(0, 15))
read2 = pd.DataFrame(df1['read2'].str.slice(0, 15))

# Join R1 and R2 together
df2 = read1.join(read2)

# Import map of barcodes to proteins
barcodes = pd.read_csv('./barcodes_run22.csv')

# Save output as df2
df2.to_csv('df2.csv')

#### PART 3: Match the extracted barcodes from the FASTQ files with the barcodes map
# Merge with read 1
# Add suffix '_x' to column names
barcodes_x = barcodes.add_suffix('_x')
df3 = pd.merge(left = df2, right = barcodes_x, left_on = "read1", right_on = "Barcode_x", how = "inner")

# Merge with read 2
# Add suffix '_y' to column names
barcodes_y = barcodes.add_suffix('_y')
df3 = pd.merge(left = df3, right = barcodes_y, left_on = "read2", right_on = "Barcode_y", how = "inner")

# Print number of rows were both read1 and read2 barcodes mapped
print('Total reads: ' + str(len(df3)) + '\n---')

# Save output as df3
df3.to_csv('df3.csv')

#### PART 4: Count the number of times each protein-protein pair occurs
# Group by Name_x (MATalpha yeast strains) and Name_y (MATa yeast strains)
description_count_series = df3.groupby(by = ['Name_x', 'Name_y']).size()

# Turn series into a data frame
description_count_df = description_count_series.to_frame(name = 'size').reset_index()

# Pivot table to M by N data frame containing a matrix of barcode pairs
df4 = pd.pivot_table(description_count_df, index = 'Name_y', columns = 'Name_x', values = 'size')

# Remove 'Name_x' and 'Name_y' headers
df4 = df4.rename_axis(None, axis = 0)
df4 = df4.rename_axis(None, axis = 1)

# Re-order columns
df4 = df4[["BFL1", "BCL2", "AAYS_263", "AAYS_265", "AAYS_356", "BCLB", "BCLw", "BCLXL", "AAYS_567", 
           "AAYS_571", "AAYS_354", "ANeg3", "AAYS_231", "ANeg1", "AAYS_566", "AAYS_232", "ANeg2"]]

# Replace NaN values with zeros
df4.fillna(0, inplace = True)

# Save output as df4
df4.to_csv('df4.csv')

#### PART 5: Visualize protein-protein interactions as a heatmap
# Transpose data frame
df4_t = df4.transpose()

# Save column names as a series
mat_a = pd.Series(df4_t.columns)
# Save index for yticks
mat_a_index = mat_a.index

# Save row names as a series
mat_alpha = pd.Series(df4_t.index)
# Save index for xticks
mat_alpha_index = mat_alpha.index

# Turn data frame into a numpy array of values
size_array = df4.values

# Create the heatmap
fig = plt.figure()
# Set figure size
fig.set_size_inches(8, 8)
ax = fig.add_subplot(1, 1, 1)
heatmap = plt.pcolor(size_array, cmap = 'Blues', edgecolors = 'black')

# Set x-axis labels, have labels appear in middle of cell
ax.set_xticks([float(n) + 0.5 for n in mat_alpha_index])
ax.set_xticklabels(mat_alpha, rotation = 90)

# Set y-axis labels, have labels appear in middle of cell
ax.set_yticks([float(n) + 0.5 for n in mat_a_index])
ax.set_yticklabels(mat_a)

# Make cells in the heatmap square
plt.axis('scaled')

# Add plot title
plt.title('Barcode pair counts')

# Label x and y axes
plt.xlabel('Mat-A')
plt.ylabel('Mat-Alpha')

# Add a colorbar
plt.colorbar()

# Show the figure
plt.show()

#### PART 6: Visualize protein-protein interactions as a Sankey diagram

# Group by Name_x (MATalpha yeast strains) and Name_y (MATa yeasts strains), keep source/target interactions with size of 0
description_count_series2 = df3.groupby(by = ['Name_x', 'Name_y']).size().unstack(fill_value=0).stack()

# Turn series into a dataframe
description_count_df2 = test.to_frame(name = 'size').reset_index()

# Generate array of a protein labels
label_a = test2['Name_x'].unique()

# Generate array of alpha protein labels
label_alpha = np.array(test2['Name_y'].head(19))

# Join a and alpha protein labels together in one array
labels = np.concatenate((label_a, label_alpha), axis = None)

# Give a proteins numbers 0-16
sources = np.arange(0, 17, 1, dtype = int)
sources2 = np.repeat(sources, 19)

# Give alpha proteins numbers 17-35
targets = np.arange(17, 36, 1, dtype = int)
targets2 = np.tile(target, 17)

# Get sizes of protein interactions
values = np.array(description_count_df2['size'])

# Plot Sankey Diagram
fig = go.Figure(data = [go.Sankey(
    node = dict(
        pad = 15,
        thickness = 20,
        line = dict(color = "black", width = 0.5),
        label = labels,
        color = "blue"
    ),
    link = dict(
      source = sources2, 
      target = target2,
      value = values
))])

fig.update_layout(
    autosize=False,
    width=500,
    height=1000
)

fig.update_layout(title_text="Barcode pair counts", font_size=10)
fig.show()
