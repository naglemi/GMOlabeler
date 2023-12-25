import os
import argparse
from PIL import ImageFont
from PIL import ImageDraw
from PIL import Image
import math
import pandas as pd
import numpy as np
import statistics
from datetime import datetime

def count_pixels(segment, rgb_from_hex):
    segment_pixel_map = segment.load()
    pixel_counter = 0
    for i in range(0, segment.size[0]):
        for j in range(0, segment.size[1]):
            if segment_pixel_map[i,j] == rgb_from_hex:
                pixel_counter = pixel_counter + 1
    return pixel_counter

def census_pixels(segment_pixel_map, segment, pixel_demographics):
    # Usage: census_pixels(segment_pixel_map = segment.load(), segment = segment, pixel_demographics = pixel_demographics)
    pixel_counts = []
    for i in range(0, pixel_demographics.shape[0]):
        n_this_color = count_pixels(segment_pixel_map = segment_pixel_map,
                                    segment = segment,
                                    rgb_from_hex = hex_to_rgb(pixel_demographics.iloc[i, 1]))
        pixel_counts.append(n_this_color)
    return pixel_counts

def hex_to_rgb(hex_color):
    hex_color = hex_color.lstrip('#')
    return tuple(int(hex_color[i:i+2], 16) for i in (0, 2, 4))

def get_row_col_for_grid_item(grid_item_in, grid_type, gmol_dir, grid_file = None):
    #print('This grid type is ' + str(grid_type))
    if grid_file==None:
        if grid_type==20:
            grid_file = pd.read_csv('explant_position_key/sections20right.txt',
                                    delimiter = '\t')
        if grid_type==12:
            grid_file = pd.read_csv('explant_position_key/sections12right.txt',
                                    delimiter = '\t')

    #print('Grid databse: ')
    #print(grid_file)
    # Column name can't be same as variable name for this to work. Hence naming of variable as grid_item_in
    #print('This grid item is ' + str(grid_item_in))
    this_grid_item_position = grid_file.loc[grid_file['grid_item'] == grid_item_in]
    #print('The row of this grid item in the grid database is ' + str(this_grid_item_position))
    row = int(this_grid_item_position['row'].iloc[0])
    col = int(this_grid_item_position['col'].iloc[0])
    return [row, col]

def determine_explant_position(grid_type,
                               grid_item,
                               gmol_dir,
                               grid_file_path = None,
                               left_edge = 14,
                               right_edge = 1406,
                               bottom_edge = 1256,
                               top_edge = 226,
                               grid_file = None):
    x_edges_cropped_size = right_edge - left_edge
    #print "X edges cropped size is " + str(x_edges_cropped_size)
    y_edges_cropped_size = bottom_edge - top_edge
    #print "Y edges cropped size is " + str(y_edges_cropped_size)

    #print "Grid item is " + str(grid_item)

    # Calculate size of each grid item
    if(grid_type==20):
        n_rows = 4
        n_cols = 6
    if(grid_type==12):
        n_rows = 3
        n_cols = 4

    if(grid_type!= 12 & grid_type!= 20):
        if(grid_file_path != None):
            stop("Error: Need to provide grid file if not using standard 12 or 20 section grid")
            # Finish up this functionality later
            grid_file = pd.read_csv(grid_file_path)
            n_rows = grid_file['row'].max
            n_cols = grid_file['col'].max

    grid_item_width = x_edges_cropped_size / n_cols
    grid_item_height = y_edges_cropped_size / n_rows

    #print "Grid item width is " + str(grid_item_width)
    #print "Grid item height is " + str(grid_item_height)

    # Find row and column of desired grid item
    row_col = get_row_col_for_grid_item(grid_item_in = grid_item,
                                        grid_type = grid_type,
                                        gmol_dir = gmol_dir,
                                        grid_file = grid_file)

    row = row_col[0]
    col = row_col[1]

    #print "Grid item " + str(grid_item) + " is at row " + str(row) + " and col " + str(col)

    # Crop to desired grid item

    # The top edge of explant is equal to the overall top edge plus the row-1 times the size of row
    cropped_top_edge = int(top_edge + ((row-1)* grid_item_height))
    # The bottom edge of explant is equal to the top edge plus the
    cropped_bottom_edge = int(top_edge + ((row)* grid_item_height))
    # The left edge of explant is equal to the overall left edge plus the col-1 times size of col
    cropped_left_edge = int(left_edge + ((col-1)* grid_item_width))
    cropped_right_edge = int(left_edge + ((col)* grid_item_width))

    #bottom_edge = (bottom_edge - ((row-1)*x_grid_item_size))
    #top_edge = (bottom_edge - (row*x_grid_item_size))
    #left_edge = (left_edge + ((col-1)*y_grid_item_size))
    #right_edge = (left_edge + (col*y_grid_item_size))

    return [cropped_left_edge,
            cropped_top_edge,
            cropped_right_edge,
            cropped_bottom_edge]

    #return [cropped_bottom_edge,
    #        cropped_top_edge,
    #        cropped_left_edge,
    #        cropped_right_edge]

def load_orient_image(image):
    image = Image.open(image).convert("RGBA")
    image = image.rotate(270, expand=True).transpose(Image.FLIP_LEFT_RIGHT)
    return(image)

def crop_to_grid_borders(object_to_crop, mode='image', grid_borders=None, verbose=False):
    if grid_borders is None:
        raise ValueError("Grid borders must be provided")

    left_edge, right_edge, bottom_edge, top_edge = grid_borders

    if verbose:
        print(f'Cropping to outer grid borders: {grid_borders}')

    if mode == 'image':
        image_cropped = object_to_crop.copy()
        image_cropped = image_cropped.crop((left_edge, top_edge, right_edge, bottom_edge))
        return image_cropped


def crop_to_explant(object_to_crop, grid_item, grid_type, gmol_dir, mode = 'image', verbose = False, grid_borders = None, grid_file = None):
        # #for i in range(1,13): # Run this to test cropping
        #    object = crop_to_explant(object_to_crop = rgb, mode = "image", grid_item = i, grid_type=12)
        #    display(object)

    if 'grid_borders' not in locals() and 'grid_borders' not in globals():

        if grid_type == 20:
            left_edge = 14
            right_edge = 1406
            bottom_edge = 1256
            top_edge = 226

        if grid_type == 12:
            left_edge = 92
            top_edge = 275
            right_edge = 1262
            bottom_edge = 1200
# "309 1251 1284 129"
    else:
        left_edge = int(grid_borders[0])
        right_edge = int(grid_borders[1])
        bottom_edge = int(grid_borders[2])
        top_edge = int(grid_borders[3])

    explant_coordinates = determine_explant_position(
    grid_type = grid_type,
    grid_item = grid_item,
    left_edge = left_edge,
    top_edge = top_edge,
    right_edge = right_edge,
    bottom_edge = bottom_edge,
    gmol_dir = gmol_dir,
    grid_file = grid_file)

    if verbose == True:
        print('Cropping to ', explant_coordinates)
    if mode == 'image':

        image_cropped = object_to_crop.copy()

        image_cropped = image_cropped.crop((explant_coordinates[0], explant_coordinates[1], explant_coordinates[2], explant_coordinates[3]))

        return(image_cropped)

    if mode == 'CLS':

        #print('Segment cropped dim is ' + str(image_cropped.size))
        #segment_array = np.array(segment_cropped)[:,:,0:3]
        #print('Segment array dim is ' + str(segment_array.shape))
        #MSE_cropped = np.sum(np.square(np.subtract(segment_array, tissue_color_tuple[0:3])), axis=2).transpose((1, 0))
        object_to_crop = object_to_crop[:, ::-1]
        CLS_cropped = object_to_crop[explant_coordinates[0]:explant_coordinates[2], explant_coordinates[1]:explant_coordinates[3]]
        #CLS_cropped = object_to_crop[explant_coordinates[1]:explant_coordinates[3], explant_coordinates[0]:explant_coordinates[2]]
        object_to_crop = object_to_crop[:, ::-1]
        CLS_cropped = CLS_cropped[:, ::-1]
        return(CLS_cropped)


def load_CLS_layer(CLS_path,  layer, format_type = 'csv'):
    if format_type == 'csv':
        CLS_data = pd.read_csv(CLS_path)
    if format_type == 'hdf':
        CLS_data = pd.read_csv(CLS_path, key = 'weights')
    CLS_data = CLS_data.sort_values(by=['rows', 'cols'])
    CLS_data = np.asarray(CLS_data[layer]).reshape((CLS_data['rows'].max(),
                                            CLS_data['cols'].max()))
    CLS_data = CLS_data[:, ::-1]
    return(CLS_data)

def CLS_to_image(CLS_matrix, cap = 100, tissue_hex='#FFFFFF'):
    #cap = 100
    #print("debugging CLS_to_image")
    # Make sure no values are below zero.
    CLS_matrix[CLS_matrix < 0] = 0
    
#     # Debugging: Check if any values are above zero (image should not be all black)
#     if np.all(CLS_matrix == 0):
#         print("Warning: All CLS_matrix values are zero after clipping negative values.")
    
    # Normalize CLS_matrix values and cap them
    CLS_matrix = np.interp(CLS_matrix, (0, cap), (0, 255)).astype(np.uint8)
    
#     # Debugging: Check if any values are above zero after normalization
#     if np.all(CLS_matrix == 0):
#         print(f"Warning: All CLS_matrix values are zero after normalization. Cap used: {cap}")
    
    # Expand the matrix to have three dimensions
    CLS_matrix_expanded = np.expand_dims(CLS_matrix, axis=2)
    
    # Apply the RGB color to the image, using the existing hex_to_rgb function
    rgb_color = hex_to_rgb(tissue_hex)  # Exclude the alpha component
    rgb_color_float = np.array(rgb_color, dtype=np.float32) / 255
    
#     # Debugging: Check that rgb_color has non-zero values
#     if np.all(np.array(rgb_color) == 0):
#         print(f"Warning: RGB color is {rgb_color}, which will result in a black image.")
    
    CLS_matrix_expanded_filled = CLS_matrix_expanded * rgb_color_float
    
    # Rotate and convert to uint8
    CLS_matrix_expanded_filled = np.rot90(CLS_matrix_expanded_filled).astype(np.uint8)
    
#     # Debugging: Check if the image data contains any non-zero values before creating the image
#     if np.all(CLS_matrix_expanded_filled == 0):
#         print("Warning: The image data is all zeros before creating the PIL Image.")
    
    img = Image.fromarray(CLS_matrix_expanded_filled, 'RGB')
    
    return img

def layers_to_image(CLS_data_layer1, CLS_data_layer2, cap1, cap2):
    #print("debugging layers_to_image")
    #CLS_data_layer1[CLS_data_layer1 < 7] = 0
    CLS_data_layer1 = np.interp(CLS_data_layer1,
                                (CLS_data_layer1.min(), cap1),
                                (0, 255))
    CLS_data_layer1 = CLS_data_layer1.astype(np.uint8)
    #print(f"Layer 1 max value after scaling: {CLS_data_layer1.max()}")

    #CLS_data_layer2[CLS_data_layer2 < 7] = 0
    CLS_data_layer2 = np.interp(CLS_data_layer2,
                                (CLS_data_layer2.min(), cap2),
                                (0, 255))
    CLS_data_layer2 = CLS_data_layer2.astype(np.uint8)
    #print(f"Layer 2 max value after scaling: {CLS_data_layer2.max()}")

    CLS_data_layer1 = np.rot90(np.expand_dims(CLS_data_layer1, axis=2))
    CLS_data_layer2 = np.rot90(np.expand_dims(CLS_data_layer2, axis=2))
    blue = np.zeros_like(CLS_data_layer1)

    CLS_matrix_expanded_filled = np.concatenate((CLS_data_layer1, CLS_data_layer2, blue), axis=2)
    img = Image.fromarray(CLS_matrix_expanded_filled, 'RGB')

    return img

def load_CLS_layers(CLS_path, layer1, layer2, format_type = 'csv'):
    # Example usage: CLS_data_layer1, CLS_data_layer2 = load_CLS_layers(CLS_path = samples['CLS_data'][15],
    #                                               layer1 = 'Chl',
    #                                               layer2 = 'DsRed')
    print('Loading CLS data from path' + str(CLS_path))
    if format_type == 'csv':
        CLS_data = pd.read_csv(CLS_path)
    if format_type == 'hdf':
        CLS_data = pd.read_hdf(CLS_path, key = 'weights')
    CLS_data = CLS_data.sort_values(by=['rows', 'cols'])
    CLS_data_layer1 = np.asarray(CLS_data[layer1]).reshape((int(CLS_data['rows'].max()),
                                            int(CLS_data['cols'].max())))
    CLS_data_layer2 = np.asarray(CLS_data[layer2]).reshape((int(CLS_data['rows'].max()),
                                            int(CLS_data['cols'].max())))
    CLS_data_layer1 = CLS_data_layer1[:, ::-1]
    CLS_data_layer2 = CLS_data_layer2[:, ::-1]
    return(CLS_data_layer1, CLS_data_layer2)

def crop_CLS_matrix(object_to_crop,
                    left_edge,#92,#14,
                    top_edge,#275,
                    right_edge,# = rgb.size[1],
                    bottom_edge):# = rgb.size[0]):

    CLS_cropped = object_to_crop[top_edge:bottom_edge, left_edge:right_edge]

    return(CLS_cropped)

def filter_CLS_desired_tissues(CLS_data_in, segment, tissue_hex, tolerance):

    #CLS_data = CLS_data_in.copy()

    segment_copy = segment.copy()
    segment_array = np.array(segment_copy)[:,:,0:3]

    tissue_color_tuple = hex_to_rgb(tissue_hex)

    CLS_cropped = crop_CLS_matrix(CLS_data_in,
                                  left_edge = 0,#92,#14,
                                  top_edge = 0,#275,
                                  right_edge = segment.size[1],
                                  bottom_edge = segment.size[0]).copy()
    MSE_cropped = np.sum(np.square(np.subtract(segment_array, tissue_color_tuple[0:3])), axis=2).transpose((1, 0))
    MSE_cropped = MSE_cropped[:, ::-1]
    CLS_cropped[MSE_cropped > tolerance] = 0
    total_tissue_pixels = len(MSE_cropped[MSE_cropped < tolerance])
    return(CLS_cropped, total_tissue_pixels)

def initialize_output_df():
    output_df = pd.DataFrame(columns = ['filename',
                                         'grid_item',
                                         'segment_hex',
                                         'intensity_threshold',
                                         'n_pixels_passing_threshold',
                                         'mean_signal',
                                         'max_signal',
                                         'total_signal',
                                         'total_pixels'],
                        index = ['observation'])
    return(output_df)

def calculate_pixels_passing_threshold(CLS_data, threshold):
    n_pixels_passing = len(CLS_data[CLS_data>threshold])
    return(n_pixels_passing)

def calculate_mean_signal(CLS_data):
    # Need to omit 0's because we zero out non-tissues of interest and only want statistics for tissues of interest
    if CLS_data[CLS_data!=0].size == 0:
        print("Warning! There is no data in this CLS matrix.")
        mean = 0
    if CLS_data[CLS_data!=0].size != 0:
        mean = CLS_data[CLS_data!=0].mean()
    return(mean)
def calculate_max_signal(CLS_data):
    # Need to omit 0's because we zero out non-tissues of interest and only want statistics for tissues of interest
    if CLS_data[CLS_data!=0].size == 0:
        print("Warning! There is no data in this CLS matrix.")
        maximum = 0
    if CLS_data[CLS_data!=0].size != 0:
        maximum = CLS_data[CLS_data!=0].max()
    return(maximum)
def calculate_total_signal(CLS_data):
    # Need to omit 0's because we zero out non-tissues of interest and only want statistics for tissues of interest
    if CLS_data[CLS_data!=0].size == 0:
        print("Warning! There is no data in this CLS matrix.")
        sum = 0
    if CLS_data[CLS_data!=0].size != 0:
        sum = CLS_data[CLS_data!=0].sum()
    return(sum)

def segment_matrix(CLS_object,
                   segment,
                   segment_filename,
                   grid_item,
                   tissue_hex,
                   threshold,
                   tolerance = 100):
    output_line = initialize_output_df()
    CLS_filtered, total_tissue_pixels = filter_CLS_desired_tissues(CLS_data_in = CLS_object,
                                                               segment = segment,
                                                               tissue_hex = tissue_hex,
                                                               tolerance = tolerance)
    output_line['filename'] = segment_filename
    output_line['grid_item'] = grid_item
    output_line['segment_hex'] = tissue_hex
    output_line['intensity_threshold'] = threshold
    #print('Total pixels ' + str(total_tissue_pixels) + ' for segment ' + tissue_hex)
    
    # Starting 5.23.21, because of a rare case when ~2 pixels are in a class and neither have signal, producing "Error!" No CLS error
    #  We will increase this past 1, up to 3
    
    if total_tissue_pixels>3:
        output_line['n_pixels_passing_threshold'] = calculate_pixels_passing_threshold(CLS_data = CLS_filtered, threshold = threshold)
        output_line['mean_signal'] = calculate_mean_signal(CLS_filtered)
        output_line['max_signal'] = calculate_max_signal(CLS_filtered)
        output_line['total_signal'] = calculate_total_signal(CLS_filtered)

    if total_tissue_pixels<=3:
        output_line['n_pixels_passing_threshold'] = 'NA'
        output_line['mean_signal'] = 'NA'
        output_line['max_signal'] = 'NA'
        output_line['total_signal'] = 'NA'

    output_line['total_pixels'] = total_tissue_pixels

    return(output_line, CLS_filtered)


def get_concat_h_multi(im_list):
    # Get the maximum height of all images
    max_height = max(im.height for im in im_list)
    
    # Total width is the sum of all image widths
    total_width = sum(im.width for im in im_list)
    
    # Create a new image with the total width and the max height
    dst = Image.new('RGB', (total_width, max_height))
    
    # Paste each image into the new image
    current_width = 0
    for im in im_list:
        dst.paste(im, (current_width, 0))
        current_width += im.width
        
    return dst

def get_concat_v(im_list):
    # Get the maximum width of all images
    max_width = max(im.width for im in im_list)
    
    # Total height is the sum of all image heights
    total_height = sum(im.height for im in im_list)
    
    # Create a new image with the max width and the total height
    dst = Image.new('RGB', (max_width, total_height))
    
    # Paste each image into the new image
    current_height = 0
    for im in im_list:
        if im.mode != 'RGB':
            print(f"Warning: Image mode is {im.mode}, converting to RGB.")
            im = im.convert('RGB')
        if not im.size[0] > 0 and im.size[1] > 0:
            print(f"Warning: Image has invalid size {im.size}.")
            continue  # Skip this image or handle appropriately
        dst.paste(im, (0, current_height))
        current_height += im.height
        
    return dst

def main(sample_df_path, threshold, layer, grid_type, gmol_dir, format_type = 'csv', grid_path = None, grid_borders = None, pixel_demographics_path = None, grid_file = None):
    
    if grid_path and os.path.exists(grid_path):
        grid = load_orient_image(grid_path)
    else:
        print("Warning: No grid image found.")

    if pixel_demographics_path:
        # Load pixel_demographics from the provided path
        pixel_demographics = pd.read_csv(pixel_demographics_path)
    else:
        # Use default pixel_demographics
        pixel_demographics = pd.DataFrame(list(zip(['Shoot', 'Callus', 'Stem', 'Background'],
                                                   ['008000', '000080', '800000', '000000'],
                                                   ['green', 'blue', 'red', 'black'])),
                                          columns=['Tissue', 'hex_code', 'color'])

    job_id = os.path.dirname(sample_df_path.split("/", 4)[4])

    target_directory = gmol_dir + '/' + job_id + '/' + layer + '/'
    print("Target directory: " + target_directory)
    
    os.makedirs(target_directory, exist_ok=True)
    #os.chdir(target_directory)

    samples = pd.read_csv(sample_df_path)

    output_df = initialize_output_df()
    print('Grid type ' + str(grid_type))

    for plate in range(0, samples.shape[0]):
        print('Loading plate ' + str(plate + 1) + ' of ' + str(samples.shape[0]))
        segment_filename = samples['segment'][plate]
        print(segment_filename)
        segment = load_orient_image(segment_filename)
        #print('Plate segment loaded')
        
        # Initialize to None
        rgb = None
        rgb_gridded = None

        # Check if 'rgb' column exists and contains a valid path
        if 'rgb' in samples.columns and pd.notna(samples['rgb'][plate]):
            rgb = load_orient_image(samples['rgb'][plate])
            #print('Plate RGB loaded')
            # Check that rgb and segment sizes are the same
            assert rgb is not None and rgb.size == segment.size, \
                    f"The dimensions of the rgb and segment images must be the same. RGB size: {rgb.size}, Segment size: {segment.size}"
            # Assert that the heights are the same
            assert rgb.size[1] == grid.size[1], "RGB and Grid images must have the same height."

            # Determine the new width based on the wider image
            new_width = max(rgb.size[0], grid.size[0])

            # Pad the narrower image
            if rgb.size[0] < new_width:
                padded_rgb = Image.new('RGB', (new_width, rgb.size[1]), (0, 0, 0))
                padded_rgb.paste(rgb, (0, 0))
                rgb_to_blend = padded_rgb
            else:
                rgb_to_blend = rgb

            if grid.size[0] < new_width:
                padded_grid = Image.new('RGB', (new_width, grid.size[1]), (0, 0, 0))
                padded_grid.paste(grid, (0, 0))
                grid_to_blend = padded_grid
            else:
                grid_to_blend = grid
           
            # Ensure that both images are in the same mode
            rgb_to_blend = rgb_to_blend.convert('RGB')
            grid_to_blend = grid_to_blend.convert('RGB')

            # Assert that rgb_to_blend and grid_to_blend have the same dimensions
            assert rgb_to_blend.size == grid_to_blend.size, \
                    f"Images to blend must have the same dimensions. RGB to blend size: {rgb_to_blend.size}, Grid to blend size: {grid_to_blend.size}"

            # Continue with operations like blending
            rgb_gridded = Image.blend(rgb_to_blend, grid_to_blend, alpha=0.5)
            
            if plate == 0:
                rgb_gridded.save(target_directory + "/rgb_and_grid_test.png")

                rgb_gridded_cropped = crop_to_grid_borders(rgb_gridded, grid_borders=grid_borders, verbose=True)

                # Save the cropped blended image
                rgb_gridded_cropped.save(target_directory + "/rgb_gridded_cropped_test.png")

        else:
            print('RGB column missing or invalid. rgb and rgb_gridded set to None.')
        
        #print(plate)
        CLS_data_layer1, CLS_data_layer2 = load_CLS_layers(CLS_path = samples['CLS_data'][plate],
        layer1 = 'Chl',
        layer2 = layer,
        format_type = format_type)

        # Patch to deal with changes in laser cutoff
        # If numpy array is less wide than RGB, we will pad it with zeros to be as wide.
        if CLS_data_layer1.shape[1] < segment.size[1]:
            print("CLS image less wide than RGB. Padding width.")
            CLS_data_layer1 = np.pad(CLS_data_layer1,
                                     [(0,
                                       0),
                                     (0,
                                      segment.size[1] - CLS_data_layer1.shape[1])],
                                    'constant')

            CLS_data_layer2 = np.pad(CLS_data_layer2,
                                     [(0,
                                      0),
                                      (0,
                                       segment.size[1] - CLS_data_layer2.shape[1])],
                                      'constant')


        # If numpy array is less long than RGB, we will pad it with zeros to be as long.
        if CLS_data_layer1.shape[0] < segment.size[0]:
            print("CLS image less long than RGB. Padding width.")
            CLS_data_layer1 = np.pad(CLS_data_layer1,
                                     [(0,
                                       segment.size[0] - CLS_data_layer1.shape[0]),
                                      (0,
                                       0)],
                                      'constant')

            CLS_data_layer2 = np.pad(CLS_data_layer2,
                                     [(0,
                                       segment.size[0] - CLS_data_layer2.shape[0]),
                                      (0,
                                       0)],
                                      'constant')

            if rgb:
                rgb = rgb.crop((0, # left
                                segment.size[1] - CLS_data_layer1.shape[1], # top
                                CLS_data_layer1.shape[0], # right
                                segment.size[1])) #bottom
                rgb_gridded = rgb_gridded.crop((0, # left
                                rgb_gridded.size[1] - CLS_data_layer1.shape[1], # top
                                CLS_data_layer1.shape[0], # right
                                rgb_gridded.size[1])) #bottom

            segment = segment.crop((0, # left
                                    segment.size[1] - CLS_data_layer1.shape[1], # top
                                    CLS_data_layer1.shape[0], # right
                                    segment.size[1])) #bottom 



        for grid_item in range(1,grid_type+1):
            #print(grid_item)
            ##################################################
            ################ CROP EVERYTHING #################
            ##################################################

            if rgb:
                rgb_cropped = crop_to_explant(object_to_crop = rgb,
                                              grid_item = grid_item,
                                              grid_type = grid_type,
                                              mode = 'image',
                                              gmol_dir = gmol_dir,
                                              grid_borders = grid_borders,
                                              grid_file = grid_file)
                rgb_gridded_cropped = crop_to_explant(object_to_crop = rgb_gridded,
                                                      grid_item = grid_item,
                                                      grid_type = grid_type,
                                                      mode = 'image',
                                                      gmol_dir = gmol_dir,
                                                      grid_borders = grid_borders,
                                                      grid_file = grid_file)
            segment_cropped = crop_to_explant(object_to_crop = segment,
                                              grid_item = grid_item,
                                              grid_type = grid_type,
                                              mode = 'image',
                                              gmol_dir = gmol_dir,
                                              grid_borders = grid_borders,
                                              grid_file = grid_file)
            CLS_cropped1 = crop_to_explant(object_to_crop = CLS_data_layer1,
                                           grid_item = grid_item,
                                           grid_type = grid_type,
                                           mode = 'CLS',
                                           gmol_dir = gmol_dir,
                                           grid_borders = grid_borders,
                                           grid_file = grid_file)
            CLS_cropped2 = crop_to_explant(object_to_crop = CLS_data_layer2,
                                           grid_item = grid_item,
                                           grid_type = grid_type,
                                           mode = 'CLS',
                                           gmol_dir = gmol_dir,
                                           grid_borders = grid_borders,
                                           grid_file = grid_file)

            ##################################################
            ####### PREPARE SOME HYPERSPECTRAL IMAGES ########
            ##################################################

            pseudo_fluor = layers_to_image(CLS_data_layer1 = CLS_cropped1,#CLS_data_layer1,
                                           CLS_data_layer2 = CLS_cropped2,#CLS_data_layer2,
                                           cap1 = 200,
                                           cap2 = 400)
            bw_fluor = CLS_to_image(CLS_matrix = CLS_cropped2,
                                    cap = 400)

            ##################################################
            ######### RUN OVER EACH TISSUE SEGMENT ###########
            ##################################################

            segment_dictionary = {}
            for i in range(0, pixel_demographics.shape[0]):
                tissue = pixel_demographics['Tissue'][i]
                #print('Running for tissue ' + tissue)
                tissue_hex = pixel_demographics['hex_code'][i]
                output_line, CLS_filtered = segment_matrix(CLS_object = CLS_cropped2,
                                                           segment = segment_cropped,
                                                           segment_filename = segment_filename,
                                                           grid_item = grid_item,
                                                           tissue_hex = tissue_hex,
                                                           threshold = threshold)
                segment_signal = CLS_to_image(CLS_matrix = CLS_filtered,
                                              cap = 400,
                                              tissue_hex = tissue_hex)
                segment_dictionary[tissue] = segment_signal
                # Check the DataFrame creation from output_line
                #print(output_line)

                # Concatenate the new row DataFrame with output_df
                #output_df = pd.concat([output_df, pd.DataFrame.from_records([output_line])])
                output_df = pd.concat([output_df, output_line], ignore_index=True)


            #print("Head of output df")
            #print(output_df.head())
            ##################################################
            ###### PREPARE IMAGE OUTPUT FOR INSPECTION #######
            ##################################################
            
            def arrange_panels(panels, max_columns=4):
                """
                Arrange a list of image panels into rows with a maximum number of columns.
                """
                rows = []
                for i in range(0, len(panels), max_columns):
                    row_panels = panels[i:i + max_columns]
                    # Debugging: Check if the row panels are not empty and have content
                    #for idx, panel in enumerate(row_panels):
#                         if panel.size == 0:
#                             print(f"Warning: Panel at index {i + idx} is empty.")
#                         else:
#                             non_black_pixels = np.sum(np.array(panel) > 0)
#                             if non_black_pixels == 0:
#                                 print(f"Warning: Panel at index {i + idx} is entirely black.")

                    rows.append(get_concat_h_multi(row_panels))
                return rows

            # Define a debugging function to check for black or empty images
            def check_image_black(img, name):
                if img.size == 0:
                    print(f"Warning: The image for {name} is empty.")
                else:
                    non_black_pixels = np.sum(np.array(img) > 0)
                    if non_black_pixels == 0:
                        print(f"Warning: The image for {name} is entirely black.")
                        
            def add_label_to_image(image, label, font_path="arial.ttf", font_size=30):
                """
                Adds a label to the top of an image.

                :param image: The PIL Image object to add the label to.
                :param label: The text label to add to the image.
                :param font_path: The path to the .ttf font file.
                :param font_size: The size of the font.
                :return: An Image object with the label.
                """
                # Create an ImageDraw object
                draw = ImageDraw.Draw(image)

                # Define the font
                try:
                    font = ImageFont.truetype(font_path, font_size)
                except IOError:
                    font = ImageFont.load_default()

                # Calculate the bounding box of the text
                text_bbox = draw.textbbox((0, 0), label, font=font)
                text_width, text_height = text_bbox[2], text_bbox[3]

                # Create a new image with space for the label
                new_image_height = image.height + text_height + 10  # 10 pixels for padding
                new_image = Image.new('RGB', (image.width, new_image_height), "white")
                new_image.paste(image, (0, text_height + 10))

                # Add the text to the new image
                text_position = ((new_image.width - text_width) // 2, 5)  # 5 pixels from the top
                draw = ImageDraw.Draw(new_image)
                draw.text(text_position, label, font=font, fill="black")

                return new_image


            # Initialize list to hold all labeled panels
            all_panels_with_labels = []

            # If rgb images are present, label and add them to the all_panels_with_labels list
            if rgb:
                labeled_rgb_cropped = add_label_to_image(rgb_cropped, "rgb_cropped")
                #check_image_black(rgb_cropped, "rgb_cropped")
                all_panels_with_labels.append(labeled_rgb_cropped)

                labeled_rgb_gridded_cropped = add_label_to_image(rgb_gridded_cropped, "rgb_gridded_cropped")
                #check_image_black(rgb_gridded_cropped, "rgb_gridded_cropped")
                all_panels_with_labels.append(labeled_rgb_gridded_cropped)

            # Label and add segment and pseudo_fluor to the all_panels_with_labels list
            #check_image_black(segment_cropped, "segment_cropped")
            labeled_segment_cropped = add_label_to_image(segment_cropped, "segment_cropped")
            all_panels_with_labels.append(labeled_segment_cropped)

            #check_image_black(pseudo_fluor, "pseudo_fluor")
            labeled_pseudo_fluor = add_label_to_image(pseudo_fluor, "pseudo_fluor")
            all_panels_with_labels.append(labeled_pseudo_fluor)

            # Label and add reporter_channel_combined to the all_panels_with_labels list
            reporter_channel_combined = bw_fluor.copy()
            #check_image_black(reporter_channel_combined, "reporter_channel_combined")
            labeled_reporter_channel_combined = add_label_to_image(reporter_channel_combined, "reporter_channel_combined")
            all_panels_with_labels.append(labeled_reporter_channel_combined)

            # Label and add all other segment images to the all_panels_with_labels list
            for key in segment_dictionary:
                #check_image_black(segment_dictionary[key], key)
                labeled_segment = add_label_to_image(segment_dictionary[key], key)
                all_panels_with_labels.append(labeled_segment)

            # Dynamically create rows with labeled images
            rows = arrange_panels(all_panels_with_labels)

            # Combine all rows vertically to get the final image
            final_image = get_concat_v(rows)

            # Save the final image
            #print("Name of segment: " + segment_filename)
            #print("Current working directory where we save logic plots: " + os.getcwd())
            output_name = segment_filename.split("/")[-1].replace('_segment_uncropped_processed.png', ('_gridspot' + str(grid_item) + '.png'))
            final_image.convert("RGBA").save(target_directory + "/" + output_name)
        
        #print("Current working directory where we save stats: " + os.getcwd())
        output_df.to_csv(target_directory + "/" + 'stats.csv')

if __name__== "__main__":
    parser = argparse.ArgumentParser()

    # Define positional arguments with nargs='?' to allow for optional named arguments
    parser.add_argument('sample_df_path', help='Path to sample dataframe')
    parser.add_argument('grid_path', help='Path to grid image')
    parser.add_argument('threshold', type=float, help='Threshold value')
    parser.add_argument('layer', help='Layer name')
    parser.add_argument('grid_type', type=int, help='Grid type')
    parser.add_argument('format_type', help='File format')
    parser.add_argument('gmol_dir', help='GMOL directory path')
    parser.add_argument('grid_borders', help='Grid borders', nargs='?')

    # Define optional named argument for the segmentation model key
    parser.add_argument('--segmentation_model_key', help='Segmentation model key', default=None)
    
    # Define optional argument for user-provided grid key (see explant position key examples on repo)
    parser.add_argument('--grid_file',
                        help='Optional argument for user-provided grid key (see explant position key examples on repo)',
                        default=None)

    # Parse the arguments
    args = parser.parse_args()
    
    # We need to parse this string to a list of integers
    if args.grid_borders:
        grid_borders = [int(value) for value in args.grid_borders.split()]

    # Call the main function with the parsed arguments
    main(
        sample_df_path=args.sample_df_path,
        grid_path=args.grid_path,
        threshold=args.threshold,
        layer=args.layer,
        grid_type=args.grid_type,
        format_type=args.format_type,
        gmol_dir=args.gmol_dir,
        grid_borders=grid_borders,
        pixel_demographics_path=args.segmentation_model_key,
        grid_file=args.grid_file
    )
