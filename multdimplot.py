#! /usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
from scipy.interpolate import griddata
from matplotlib.colors import LinearSegmentedColormap
import os
import glob

print("RAPTOR multi-frequency imaging visualization script \n")

# Get a list of all image data files
data_files = sorted(glob.glob('output/img_data_*_*_60.00.dat'))

if not data_files:
    print("No data files found!")
    exit(1)

# Extract the frequencies from the filenames
frequencies = []
for file in data_files:
    parts = file.split('_')
    freq = float(parts[3])
    frequencies.append(freq)

print(f"Found {len(frequencies)} frequency files: {frequencies}")

# Create output directory if it doesn't exist
os.makedirs('multi_freq_images', exist_ok=True)

# Define a color for each frequency on a rainbow spectrum
colors = plt.cm.rainbow(np.linspace(0, 1, len(frequencies)))

# Store all processed images for the merged visualization
all_images = []
img_size = 0

# Process each frequency file
for i, filename in enumerate(data_files):
    print(f"\nProcessing {filename}...")

    # Read header
    with open(filename) as f:
        line = f.readline()
        header = line.split()

    print("IMAGE INFORMATION")
    print("Image size", int(header[0]), "x", int(header[1]))
    print("Flux of", float(header[2]), "at frequency", float(header[3]))

    img_size = int(header[0])
    img_size2 = int(header[1])

    # Read data
    data = np.loadtxt(filename, skiprows=1) + 1e-20

    # Process data
    if(img_size != img_size2):
        xi = np.linspace(0, img_size, img_size2)
        yi = np.linspace(0, img_size, img_size2)
        xi = np.tile(xi, img_size2)
        yi = np.repeat(yi, img_size2)
        x = np.linspace(0, img_size, img_size)
        y = np.linspace(0, img_size, img_size)
        x = np.tile(x, img_size)
        y = np.repeat(y, img_size)
        img = griddata((x, y), data, (xi, yi), method='cubic', rescale=True)
    else:
        img = data

    img = np.reshape(data, (-1, img_size))
    img = np.transpose(img)
    img = np.flipud(img)

    # Normalize image for individual plotting
    norm_img = np.sqrt(img / np.max(img))

    # Store the original image for merged visualization
    all_images.append(img)

    # Plot individual frequency image
    plt.clf()
    plt.close('all')

    fig = plt.figure(figsize=(10, 10))
    halfrange = 20
    ax = fig.add_subplot(111)

    im = plt.imshow(norm_img, interpolation='Nearest', cmap='magma',
                   extent=[-halfrange, halfrange, -halfrange, halfrange])

    # Add frequency label
    plt.title(f"Frequency: {frequencies[i]:.2e} Hz", fontsize=12)
    plt.colorbar(label='Normalized Intensity')

    # Save individual image
    output_file = f"multi_freq_images/img_freq_{frequencies[i]:.2e}.png"
    plt.savefig(output_file, bbox_inches='tight', dpi=150)
    print(f"Saved {output_file}")

# Create merged image with different colors for each frequency
plt.clf()
plt.close('all')

fig = plt.figure(figsize=(12, 10), facecolor='black')
halfrange = 20
ax = fig.add_subplot(111)
ax.set_facecolor('black')

# Create a merged image
merged_img = np.zeros((img_size, img_size, 4))  # RGBA image

# Add each frequency with its own color
for i, img in enumerate(all_images):
    # Normalize each image
    norm_img = img / np.max(img)

    # Create a colored RGBA version of this frequency's image
    color_img = np.zeros((img_size, img_size, 4))

    # Use the color from our color map for this frequency
    r, g, b, a = colors[i]

    # Fill in the RGBA values
    for c in range(3):  # RGB channels
        color_img[:, :, c] = norm_img * colors[i, c]

    # Set alpha based on intensity
    color_img[:, :, 3] = np.sqrt(norm_img) * 0.8  # Use sqrt for better visibility

    # Add to the merged image (simple addition)
    merged_img += color_img

# Normalize the merged image if needed
merged_img = np.clip(merged_img, 0, 1)

# Plot merged image
im = plt.imshow(merged_img, interpolation='Nearest',
               extent=[-halfrange, halfrange, -halfrange, halfrange])

# Create legend for frequencies
for i, freq in enumerate(frequencies):
    plt.plot([], [], color=colors[i], label=f"{freq:.2e} Hz", linewidth=3)

legend = plt.legend(loc='upper right', title="Frequencies")
legend.get_frame().set_facecolor('black')
legend.get_frame().set_edgecolor('white')
plt.setp(legend.get_texts(), color='white')
plt.setp(legend.get_title(), color='white')
plt.title("Merged Multi-Frequency Image", fontsize=14, color='white')
plt.xlabel(r"Distance ($R_{\mathrm{g}})$", color='white')
plt.ylabel(r"Distance ($R_{\mathrm{g}})$", color='white')
ax.tick_params(axis='x', colors='white')
ax.tick_params(axis='y', colors='white')

# Save merged image
merged_output_file = "multi_freq_images/merged_image.png"
plt.savefig(merged_output_file, bbox_inches='tight', dpi=200, facecolor='black')
print(f"\nSaved merged image: {merged_output_file}")

print("\nDone!")
