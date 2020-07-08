#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright: (c) 2019 Centre National d'Etudes Spatiales

"""
This module contains the main functions to generate the dataset
"""

from argparse import ArgumentParser
from json_checker import CheckerError
import logging
import os
import json
import subprocess
from osgeo import gdal


import configuration as cfg
from create_disp_epi import create_disp_epi
from validation import validation
from validation.validation_ratio_disp import validate_ratio_disp
from lidar_raster_registration import lidar_raster_registration
from cars import utils
import processing as proc
from georeference_dfc_files import georeference_images
from processing import get_epsg


def run(cfg, out_dir, verbose):
    """
    Pipeline to generate disparity maps in epipolar geometry
    :param cfg: configuration
    :param verbose: verbose mode
    :return:
    """

    # out statistics
    out_stats_filename = os.path.join(out_dir, 'output_statistics.json')
    if os.path.exists(out_stats_filename):
        with open(out_stats_filename, 'r') as f:
            out_stats = json.load(f)
    else:
        out_stats = {}
        
    
        
    # Apply lidar registration if needed
    if 'processing' in cfg:
        if "lidar_raster_registration" in cfg['processing']:
            #  Georeference images if needed
            if "dataset" in cfg['input']:
                if cfg['input']['dataset'] == 'DFC':
                    # Georeference lidar
                    ds = gdal.Open(cfg['processing']['lidar_raster_registration']["unregistered_lidar"])
                    root_ext = cfg['processing']['lidar_raster_registration']["unregistered_lidar"].split(os.extsep)
                    txt_file = root_ext[0] + '.txt'
                    if ds.GetProjection() == '':
                        out_dsm_georef = os.path.join(out_dir, 'DSM_georef.tif')
                        # Get the epsg code from the image
                        epsg = get_epsg(cfg['input']['img1'])
                        georeference_images(cfg['processing']['lidar_raster_registration']["unregistered_lidar"], txt_file,out_dsm_georef, epsg)
                        cfg['processing']['lidar_raster_registration']["unregistered_lidar"] = out_dsm_georef
                    # Georeference classification
                    if 'classif' in cfg['input']:
                        ds = gdal.Open(cfg['input']['classif'])
                        if ds.GetProjection() == '':
                            # Get the epsg code from the image
                            epsg = get_epsg(cfg['input']['img1'])
                            out_classif_georef = os.path.join(out_dir, 'CLASSIF_georef.tif')
                            georeference_images(cfg['input']['classif'], txt_file, out_classif_georef, epsg)
                            cfg['input']['classif'] = out_classif_georef
    
            # Apply lidar raster registration, and update the configuration
            cfg, out_stats['lidar_registration'] = lidar_raster_registration(cfg, out_dir, "median")
            # Images are not consistent, exit the program
            if cfg is None:
                return 0
        
    
        strategy = "superimpose"
        if "strategy" not in cfg["processing"]["create_disp_epi"]:
            logger.INFO("default strategy {}".format(strategy))
        else:
            strategy = cfg["processing"]["create_disp_epi"]["strategy"]
    
        # Generate disparity maps in epipolar geometry
        out_stats['create_disp'], epipolar_size_x, epipolar_size_y = create_disp_epi(cfg, out_dir, strategy)


    # Validate the disparity map if needed
    if 'validation' in cfg:
        if cfg['validation']['fractional_histogram']:
            validation.fractional_histogram(os.path.join(out_dir, 'left_epipolar_disp.tif'),
                                            os.path.join(out_dir, 'valid_disp.tif'), bins=20,
                                            title="Histogram of disparity fractional parts",
                                            savefig=os.path.join(out_dir, 'Validation_histogram_ref_disparity.tif'))
        if cfg['validation']['right_image_resampling']:
            out_stats['right_to_left_resampling_image_diff'] = validation.right_image_resampling(os.path.join(out_dir, 'left_epipolar_disp.tif'),
                                              os.path.join(out_dir, 'right_epipolar_image.tif'),
                                              os.path.join(out_dir, 'left_epipolar_image.tif'),
                                              os.path.join(out_dir, "valid_disp.tif"),
                                              os.path.join(out_dir, 'Validation_right_epipolar_resampled_image.tif'),
                                              os.path.join(out_dir, 'Validation_difference_left_right_resampled.tif'))

        if cfg['validation']['validate_ratio_disp']:
            out_stats['validate_ratio_disp'] = validate_ratio_disp(os.path.join(out_dir, 'left_monoband.tif'),
                                                                   os.path.join(out_dir, 'grid_left.tif'),
                                                                   os.path.join(out_dir, 'right_monoband.tif'),
                                                                   os.path.join(out_dir, 'right_epipolar_grid.tif'),
                                                                   epipolar_size_x, epipolar_size_y)

        if cfg['validation']['keypoints']['activate']:
            out_stats['validation_sift_ground'] = validation.keypoints_validation(os.path.join(out_dir, 'left_epipolar_image.tif'),
                                            os.path.join(out_dir, 'right_epipolar_image.tif'),
                                            os.path.join(out_dir, 'left_epipolar_classif.tif'),
                                            os.path.join(out_dir, 'left_epipolar_disp.tif'),
                                            os.path.join(out_dir, "valid_disp.tif"),
                                            os.path.join(out_dir, 'left_ratio_disparity_altitude.tif'),
                                            cfg['processing']['create_disp_epi']['registered_lidar'],
                                            os.path.join(out_dir, 'left_bias_disparity_altitude.tif'),
                                            'sift',
                                            10.0,
                                            2)
            out_stats['validation_sift_building'] = validation.keypoints_validation(os.path.join(out_dir, 'left_epipolar_image.tif'),
                                            os.path.join(out_dir, 'right_epipolar_image.tif'),
                                            os.path.join(out_dir, 'left_epipolar_classif.tif'),
                                            os.path.join(out_dir, 'left_epipolar_disp.tif'),
                                            os.path.join(out_dir, "valid_disp.tif"),
                                            os.path.join(out_dir, 'left_ratio_disparity_altitude.tif'),
                                            cfg['processing']['create_disp_epi']['registered_lidar'],
                                            os.path.join(out_dir, 'left_bias_disparity_altitude.tif'),
                                            'sift',
                                            10.0,
                                            6)
            out_stats['validation_surf_ground'] = validation.keypoints_validation(os.path.join(out_dir, 'left_epipolar_image.tif'),
                                            os.path.join(out_dir, 'right_epipolar_image.tif'),
                                            os.path.join(out_dir, 'left_epipolar_classif.tif'),
                                            os.path.join(out_dir, 'left_epipolar_disp.tif'),
                                            os.path.join(out_dir, "valid_disp.tif"),
                                            os.path.join(out_dir, 'left_ratio_disparity_altitude.tif'),
                                            cfg['processing']['create_disp_epi']['registered_lidar'],
                                            os.path.join(out_dir, 'left_bias_disparity_altitude.tif'),
                                            'surf',
                                            10.0,
                                            2)
            out_stats['validation_surf_building'] = validation.keypoints_validation(os.path.join(out_dir, 'left_epipolar_image.tif'),
                                            os.path.join(out_dir, 'right_epipolar_image.tif'),
                                            os.path.join(out_dir, 'left_epipolar_classif.tif'),
                                            os.path.join(out_dir, 'left_epipolar_disp.tif'),
                                            os.path.join(out_dir, "valid_disp.tif"),
                                            os.path.join(out_dir, 'left_ratio_disparity_altitude.tif'),
                                            cfg['processing']['create_disp_epi']['registered_lidar'],
                                            os.path.join(out_dir, 'left_bias_disparity_altitude.tif'),
                                            'surf',
                                            10.0,
                                            6)
    
    if 'epipolar_registration_analysis' in cfg['validation']:
        if cfg['validation']['epipolar_registration_analysis']["activate"]:
            singularity_image =  cfg['validation']['epipolar_registration_analysis']['singularity_image']
            ra_script = cfg['validation']['epipolar_registration_analysis']['ra_script']
            medicis_props = cfg['validation']['epipolar_registration_analysis']['medicis_props']
            left_image = os.path.join(out_dir, 'left_epipolar_image.tif')
            right_image = os.path.join(out_dir, 'right_epipolar_image.tif')
            out_stats['epipolar_registration_analysis'] = validation.epipolar_registration_analysis(os.path.join(out_dir, left_image),
                                              right_image,
                                              singularity_image,
                                              ra_script,
                                              os.path.join(out_dir, "temp"),
                                              medicis_props,
                                              os.path.join(out_dir, 'Validation_epipolar_resampling_analysis.png'),
                                              cfg['validation']['epipolar_registration_analysis']['resampling_factor'])
            out_stats['epipolar_registration_analysis']['left_image'] = left_image
            out_stats['epipolar_registration_analysis']['right_image'] = right_image
            out_stats['epipolar_registration_analysis']['resampling_factor'] = cfg['validation']['epipolar_registration_analysis']['resampling_factor']
            out_stats['epipolar_registration_analysis']['band'] = cfg['validation']['epipolar_registration_analysis']['band']
    with open(out_stats_filename, 'w') as f:
        json.dump(out_stats, f, indent=4)

if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument("-v", "--verbose", action="store_true", dest="verbose",
                        default=False, help="Verbose mode")
    parser.add_argument("injson", help="Input json file")
    parser.add_argument("outdir", help="Output directory")
    args = parser.parse_args()
    injson = cfg.read_json_file(args.injson)

    # Check configuration with respect to schema
    try:
        configuration = cfg.check_json(injson)
    except CheckerError as e:
        logging.warning(
            "Input json does not comply with schema: {}".format(e))
        exit()
    logger = logging.getLogger()
    if args.verbose:
        logger.setLevel(logging.INFO)
    else:
        logger.setLevel(logging.WARNING)

    run(configuration, os.path.abspath(args.outdir), args.verbose)
