import os

import sys
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

def hex_to_rgb(hex_code):
    rgbt_values = tuple(int(hex_code[i:i+2], 16) for i in (0, 2, 4)) + tuple((255,))
    return rgbt_values

def get_row_col_for_grid_item(grid_item_in, grid_type, gmol_dir, grid_file = None):
    #print('This grid type is ' + str(grid_type))
    if(grid_type==20):
        grid_file = pd.read_csv(gmol_dir + '/explant_position_key/sections20right.txt',
                               delimiter = '\t')
    if(grid_type==12):
        grid_file = pd.read_csv(gmol_dir + '/explant_position_key/sections12right.txt',
                               delimiter = '\t')

    #print('Grid databse: ')
    #print(grid_file)
    # Column name can't be same as variable name for this to work. Hence naming of variable as grid_item_in
    #print('This grid item is ' + str(grid_item_in))
    this_grid_item_position = grid_file.loc[grid_file['grid_item'] == grid_item_in]
    #print('The row of this grid item in the grid database is ' + str(this_grid_item_position))
    row = int(this_grid_item_position['row'])
    col = int(this_grid_item_position['col'])
    return [row, col]

def determine_explant_position(grid_type,
                               grid_item,
			       gmol_dir,
                               grid_file_path = None,
                               left_edge = 14,
                               right_edge = 1406,
                               bottom_edge = 1256,
                               top_edge = 226):
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
    row_col = get_row_col_for_grid_item(grid_item_in = grid_item, grid_type = grid_type, gmol_dir = gmol_dir)
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


def crop_to_explant(object_to_crop, grid_item, grid_type, gmol_dir, mode = 'image', verbose = False:
        # #for i in range(1,13): # Run this to test cropping
        #    object = crop_to_explant(object_to_crop = rgb, mode = "image", grid_item = i, grid_type=12)
        #    display(object)

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

        explant_coordinates = determine_explant_position(grid_type = grid_type, grid_item = grid_item, # Not sure why I need -1
                                                        left_edge = left_edge,#14,
                                                        top_edge = top_edge,
                                                        right_edge = right_edge,
                                                        bottom_edge = bottom_edge
                                                        gmol_dir = gmol_dir)
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


def load_CLS_layer(CLS_path,  layer, format = 'csv'):
    if format == 'csv':
        CLS_data = pd.read_csv(CLS_path)
    if format == 'hdf':
        CLS_data = pd.read_csv(CLS_path, key = 'weights')
    CLS_data = CLS_data.sort_values(by=['rows', 'cols'])
    CLS_data = np.asarray(CLS_data[layer]).reshape((CLS_data['rows'].max(),
                                            CLS_data['cols'].max()))
    CLS_data = CLS_data[:, ::-1]
    return(CLS_data)


def CLS_to_image(CLS_matrix, cap, mode = 'opaque', match_size=True, color='white'):
    # Usage:
    # image_out = CLS_to_image(CLS_matrix = load_CLS_layer(CLS_path = samples['CLS_data'][15],
    #                                                 layer = 'DsRed'),
    #                     cap = 10)
    CLS_matrix[CLS_matrix < 0] = 0
    CLS_matrix = np.interp(CLS_matrix,
                             (0,
                              #CLS_matrix.max()),
                              cap),
                             (0,
                              255)).astype(int)
    CLS_matrix_expanded = np.rot90(np.expand_dims(CLS_matrix, axis=2))
    if mode == 'opaque':
        empty_channel =  np.zeros((CLS_matrix_expanded.shape[0], CLS_matrix_expanded.shape[1], 1), dtype=np.uint8)
        if color=='white':
            CLS_matrix_expanded_filled = np.concatenate((CLS_matrix_expanded,
                                                         CLS_matrix_expanded,
                                                         CLS_matrix_expanded), axis=2).astype(np.uint8)
        if color=='red':
            CLS_matrix_expanded_filled = np.concatenate((CLS_matrix_expanded,
                                                         empty_channel,
                                                         empty_channel), axis=2).astype(np.uint8)
        if color=='green':
            CLS_matrix_expanded_filled = np.concatenate((empty_channel,
                                                         CLS_matrix_expanded,
                                                         empty_channel), axis=2).astype(np.uint8)
        if color=='blue':
            CLS_matrix_expanded_filled = np.concatenate((empty_channel,
                                                         empty_channel,
                                                         CLS_matrix_expanded), axis=2).astype(np.uint8)
        img = Image.fromarray(CLS_matrix_expanded_filled, 'RGB')
    if mode == 'transparent': # Don't think this mode is working yet
        blue_red_transparency =  np.zeros((CLS_matrix_expanded.shape[0], CLS_matrix_expanded.shape[1], 3), dtype=np.uint8)
        CLS_matrix_expanded_filled = np.concatenate((CLS_matrix_expanded, blue_red_transparency), axis=2).astype(np.uint8)
    if match_size==True:
        img = img.crop((0, 0, rgb.size[0], rgb.size[1]))

    return(img)


def load_CLS_layers(CLS_path, layer1, layer2, format = 'csv'):
    # Example usage: CLS_data_layer1, CLS_data_layer2 = load_CLS_layers(CLS_path = samples['CLS_data'][15],
    #                                               layer1 = 'Chl',
    #                                               layer2 = 'DsRed')
    print('Loading CLS data from path' + str(CLS_path))
    if format == 'csv':
        CLS_data = pd.read_csv(CLS_path)
    if format == 'hdf':
        CLS_data = pd.read_hdf(CLS_path, key = 'weights')
    CLS_data = CLS_data.sort_values(by=['rows', 'cols'])
    CLS_data_layer1 = np.asarray(CLS_data[layer1]).reshape((int(CLS_data['rows'].max()),
                                            int(CLS_data['cols'].max())))
    CLS_data_layer2 = np.asarray(CLS_data[layer2]).reshape((int(CLS_data['rows'].max()),
                                            int(CLS_data['cols'].max())))
    CLS_data_layer1 = CLS_data_layer1[:, ::-1]
    CLS_data_layer2 = CLS_data_layer2[:, ::-1]
    return(CLS_data_layer1, CLS_data_layer2)

def layers_to_image(CLS_data_layer1, CLS_data_layer2, cap1, cap2, match_size=True):
    # Usage:
    # image_out = CLS_to_image(CLS_matrix = load_CLS_layer(CLS_path = samples['CLS_data'][15],
    #                                                 layer = 'DsRed'),
    #                     cap = 10)
    CLS_data_layer1 = CLS_data_layer1.copy()
    CLS_data_layer2 = CLS_data_layer2.copy()
    CLS_data_layer1[CLS_data_layer1 < 7] = 0
    CLS_data_layer1 = np.interp(CLS_data_layer1,
                             (0,
                              #CLS_matrix.max()),
                              cap1),
                             (0,
                              255)).astype(int)
    CLS_data_layer1 = np.rot90(np.expand_dims(CLS_data_layer1, axis=2))
    blue =  np.zeros((CLS_data_layer1.shape[0], CLS_data_layer1.shape[1], 1), dtype=np.uint8)

    CLS_data_layer2[CLS_data_layer2 < 7] = 0
    CLS_data_layer2 = np.interp(CLS_data_layer2,
                             (0,
                              #CLS_matrix.max()),
                              cap2),
                             (0,
                              255)).astype(int)
    CLS_data_layer2 = np.rot90(np.expand_dims(CLS_data_layer2, axis=2))

    CLS_matrix_expanded_filled = np.concatenate((CLS_data_layer1, CLS_data_layer2, blue), axis=2).astype(np.uint8)
    img = Image.fromarray(CLS_matrix_expanded_filled, 'RGB')
    if match_size==True:
        img = img.crop((0, 0, rgb.size[0], rgb.size[1]))

    return(img)

def superimpose_grid(image, grid_path):
    # DEPRECATED
    # Usage: image_gridded = superimpose_grid(image_out, grid_type=12)
    if(grid_type==12):
        grid_path = '/scratch2/NSF_GWAS/macroPhor_Array/grids/grids_left_facing_125208_1_0_1_rgb_processed.jpg'
    if(grid_type==20):
        grid_path = '/scratch2/NSF_GWAS/macroPhor_Array/grids/grids_left_facing_125201_0_0_0_rgb_processed.jpg'
    grid = Image.open(grid_path).convert("RGB")
    #print('Size of grid and image, respectively')
    #print(grid.size)
    #print(image.size)
    grid = grid.rotate(270, expand=True).transpose(Image.FLIP_LEFT_RIGHT)
    #grid = grid.transpose(Image.FLIP_LEFT_RIGHT)
    #print(grid.size)
    output = Image.blend(image.convert("RGB"), grid, alpha=0.5)
    return(output)

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
    mean = CLS_data[CLS_data!=0].mean()
    return(mean)
def calculate_max_signal(CLS_data):
    # Need to omit 0's because we zero out non-tissues of interest and only want statistics for tissues of interest
    if CLS_data[CLS_data!=0].size == 0:
        sys.exit("Error! There is no data in this CLS matrix.")
    maximum = CLS_data[CLS_data!=0].max()
    return(maximum)
def calculate_total_signal(CLS_data):
    # Need to omit 0's because we zero out non-tissues of interest and only want statistics for tissues of interest
    sum = CLS_data[CLS_data!=0].sum()
    return(sum)

def segment_matrix(CLS_object,
                         segment,
                         segment_filename,
                         grid_item,
                         tissue_hex,
                         threshold,
                         tolerance = 25000):
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
    if total_tissue_pixels>1:
        output_line['n_pixels_passing_threshold'] = calculate_pixels_passing_threshold(CLS_data = CLS_filtered, threshold = threshold)
        output_line['mean_signal'] = calculate_mean_signal(CLS_filtered)
        output_line['max_signal'] = calculate_max_signal(CLS_filtered)
        output_line['total_signal'] = calculate_total_signal(CLS_filtered)

    if total_tissue_pixels<=1:
        output_line['n_pixels_passing_threshold'] = 'NA'
        output_line['mean_signal'] = 'NA'
        output_line['max_signal'] = 'NA'
        output_line['total_signal'] = 'NA'

    output_line['total_pixels'] = total_tissue_pixels

    return(output_line, CLS_filtered)


def get_concat_h(im1, im2): # Thanks to https://note.nkmk.me/en/python-pillow-concat-images/ for this
    dst = Image.new('RGB', (im1.width + im2.width, im1.height))
    dst.paste(im1, (0, 0))
    dst.paste(im2, (im1.width, 0))
    return dst

def get_concat_v(im1, im2): # Thanks to https://note.nkmk.me/en/python-pillow-concat-images/ for this
    dst = Image.new('RGB', (im1.width, im1.height + im2.height))
    dst.paste(im1, (0, 0))
    dst.paste(im2, (0, im1.height))
    return dst

def main(sample_df_path, grid, threshold, layer, grid_type, gmol_dir, format = 'csv'):

    pixel_demographics = pd.DataFrame(list(zip(['Shoot', 'Callus', 'Stem', 'Background'],
                                               ['00CC11', '0006CC', 'CC0000', '000000'],
                                               ['green', 'blue', 'red', 'black'])),
                                     columns = ['Tissue', 'hex_code', 'color'])

    job_id = os.path.dirname(sample_df_path.split("/", 4)[4])

    target_directory = gmol_dir + '/output/' + job_id
    os.makedirs(target_directory, exist_ok=True)
    os.chdir(target_directory)

    samples = pd.read_csv(sample_df_path)

    output_df = initialize_output_df()
    print('Grid type ' + str(grid_type))

    for plate in range(0, samples.shape[0]):
        print('Loading plate ' + str(plate) + ' of ' + str(samples.shape[0]))
        segment_filename = samples['segment'][plate]
        print(segment_filename)
        segment = load_orient_image(segment_filename)
        print('Plate segment loaded')
        rgb = load_orient_image(samples['rgb'][plate])
        print('Plate RGB loaded')
        rgb_gridded = Image.blend(rgb, grid, alpha=0.5)
        #print(plate)
        CLS_data_layer1, CLS_data_layer2 = load_CLS_layers(CLS_path = samples['CLS_data'][plate],
        layer1 = 'Chl',
        layer2 = layer,
        format = format)
        for grid_item in range(1,grid_type+1):

            ##################################################
            ################ CROP EVERYTHING #################
            ##################################################

            rgb_cropped = crop_to_explant(object_to_crop = rgb,
                                              grid_item = grid_item,
                                              grid_type = grid_type,
                                              mode = 'image',
					      gmol_dir = gmol_dir)
            rgb_gridded_cropped = crop_to_explant(object_to_crop = rgb_gridded,
                                              grid_item = grid_item,
                                              grid_type = grid_type,
                                              mode = 'image',
					      gmol_dir = gmol_dir)
            segment_cropped = crop_to_explant(object_to_crop = segment,
                                              grid_item = grid_item,
                                              grid_type = grid_type,
                                              mode = 'image',
					      gmol_dir = gmol_dir)
            CLS_cropped1 = crop_to_explant(object_to_crop = CLS_data_layer1,
                                              grid_item = grid_item,
                                              grid_type = grid_type,
                                              mode = 'CLS',
					      gmol_dir = gmol_dir)
            CLS_cropped2 = crop_to_explant(object_to_crop = CLS_data_layer2,
                                              grid_item = grid_item,
                                              grid_type = grid_type,
                                              mode = 'CLS',
					      gmol_dir = gmol_dir)

            ##################################################
            ####### PREPARE SOME HYPERSPECTRAL IMAGES ########
            ##################################################

            pseudo_fluor = layers_to_image(CLS_data_layer1 = CLS_cropped1,#CLS_data_layer1,
                                  CLS_data_layer2 = CLS_cropped2,#CLS_data_layer2,
                                  cap1 = 40,
                                  cap2 = 10,
                                  match_size = False)
            bw_fluor = CLS_to_image(CLS_matrix = CLS_cropped2,
                                     cap = 10,
                                     match_size=False)

            ##################################################
            ######### RUN OVER EACH TISSUE SEGMENT ###########
            ##################################################

            segment_dictionary = {}
            for i in range(0, pixel_demographics.shape[0]-1):
                tissue = pixel_demographics['Tissue'][i]
                #print('Running for tissue ' + tissue)
                tissue_hex = pixel_demographics['hex_code'][i]
                color = pixel_demographics['color'][i]
                output_line, CLS_filtered = segment_matrix(CLS_object = CLS_cropped2,
                                                       segment = segment_cropped,
                                                       segment_filename = segment_filename,
                                                       grid_item = grid_item,
                                                       tissue_hex = tissue_hex,
                                                       threshold = threshold)
                segment_signal = CLS_to_image(CLS_matrix = CLS_filtered,
                                     cap = 10,
                                     match_size=False,
                                     color = color)
                segment_dictionary[tissue] = segment_signal
                output_df = output_df.append(output_line)

            ##################################################
            ###### PREPARE IMAGE OUTPUT FOR INSPECTION #######
            ##################################################

            rgb_with_without_grid = get_concat_h(rgb_cropped, rgb_gridded_cropped)
            segment_and_pseudo_fluor = get_concat_h(segment_cropped, pseudo_fluor)
            top_row = get_concat_h(rgb_with_without_grid, segment_and_pseudo_fluor)
            reporter_channel_combined = bw_fluor.copy()
            for key in segment_dictionary:
                #print(key)
                reporter_channel_combined = get_concat_h(reporter_channel_combined, segment_dictionary[key])
            both_rows = get_concat_v(top_row, reporter_channel_combined)
            #display(both_rows)
            output_name = segment_filename.split("/")[-1].replace('_segment_uncropped_processed.jpg',
                                                          ('_gridspot'+str(grid_item)+'.png'))
            #print('Saving ' + str(output_name))
            both_rows.convert("RGBA").save(output_name)

        output_df.to_csv('stats.csv')

if __name__== "__main__":
    #print(sys.argv)
    #print('sys.argv[5] is ' + sys.argv[5])
    main(sample_df_path = sys.argv[1],
	 grid = load_orient_image(sys.argv[2]),
	 threshold = float(sys.argv[3]),
	 layer = sys.argv[4],
	 grid_type = int(sys.argv[5]),
	 format = str(sys.argv[6]),
	 gmol_dir = str(sys.argv[7]))
