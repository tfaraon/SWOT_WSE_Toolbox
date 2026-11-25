import numpy as np
import pandas as pd
import rasterio
from rasterio.transform import rowcol
from rasterio.warp import transform
from . import SWOT_loader as sload

import json
import os
import re
from datetime import datetime
from tqdm import tqdm

from concurrent.futures import ProcessPoolExecutor, as_completed
import multiprocessing
from functools import partial


def load_mnt(mnt_path):
    """
    Ouvre un MNT avec rasterio et retourne l'objet rasterio.io.DatasetReader.
    """
    return rasterio.open(mnt_path)

def get_mnt_value(mnt, lon, lat):
    """
    Retourne la valeur du MNT à la position (lon, lat).
    """
    
    # S'assurer que les coordonnées sont dans le bon CRS
    lon_lat = np.array([[lon, lat]])
    if mnt.crs.to_string() != 'EPSG:4326':
        # transformation depuis EPSG:4326 vers le CRS du MNT
        lon_lat = np.array(transform('EPSG:4326', mnt.crs, [lon], [lat])).T

    # Convertir coordonnées géographiques en indices ligne/colonne
    row, col = rowcol(mnt.transform, lon_lat[0][0], lon_lat[0][1])
    
    # Lire la valeur dans la bande 1
    return mnt.read(1)[row, col]

def correct_SWOT_heights(mnt, swot_dir, lon, lat): 
    """
    Corrige la hauteur SWOT par rapport à l'altitude du MNT à une coordonnée donnée.
    """
    try:
        mnt_value = get_mnt_value(mnt, lon, lat)
    except IndexError:
        print("Coordonnées en dehors des limites du MNT.")
        return None

    # Extraire la série temporelle SWOT
    wse_timeseries = extract_wse_timeseries(swot_dir, lon, lat, filter_bound=True, filter_outliers=True)
    
    if wse_timeseries.empty:
        print("Série temporelle SWOT vide à ces coordonnées.")
        return None

    wse_value = wse_timeseries['wse'].mean()
    diff = mnt_value - wse_value
    corrected_wse = wse_value + diff

    return wse_value, corrected_wse


def compute_qual_percentages(data_dict):
    """
    Calcule le pourcentage de valeurs wse_qual pour chaque clé d'un dictionnaire.
    
    Parameters
    ----------
    data_dict : dict
        Dictionnaire de séries temporelles (DataFrames) avec une colonne 'qual'
    
    Returns
    -------
    pd.DataFrame
        Tableau avec les pourcentages des qualités 0, 1, 2, 3 pour chaque clé
    """
    qual_levels = [0, 1, 2, 3]
    result = []

    for key, df in data_dict.items():
        if df is None or df.empty or 'qual' not in df.columns:
            print(f"Warning: Série vide ou invalide pour {key}")
            continue
        counts = df['qual'].value_counts(normalize=True) * 100
        row = {str(qual): counts.get(qual, 0.0) for qual in qual_levels}
        row['Label'] = key
        result.append(row)
    
    result_df = pd.DataFrame(result).set_index('Label')
    result_df = result_df[['0', '1', '2', '3']]

    

    return result_df


def extract_wse_timeseries(directory_path, lon, lat, filter_pass=None, filter_resolution=None,
                           filter_tile=None, buffer_size=0, filter_outliers=False, filter_bound=False,
                           wse_qual_filter=None, mnt_corr=None, debug=False, **station_options):
    """Notebook interface. Use station() to include rejected acquisitions."""
    return extract_wse_timeseries_parallel(directory_path, lon, lat, filter_pass, filter_resolution,
        filter_tile, buffer_size, filter_outliers, filter_bound, wse_qual_filter, mnt_corr, debug,
        n_workers=1, **station_options)


def extract_wse_timeseries_parallel(directory_path, lon, lat, filter_pass=None, filter_resolution=None,
                                    filter_tile=None, buffer_size=0, filter_outliers=False, filter_bound=False,
                                    wse_qual_filter=None, mnt_corr=None, debug=False, n_workers=None,
                                    **station_options):
    from .stations import station
    options = dict(radius=buffer_size, levels=wse_qual_filter, quality_mode="centre",
                   resolution=filter_resolution, orbit=filter_pass, tile=filter_tile,
                   bounds=(-16, 6) if filter_bound else None,
                   iqr="reject" if filter_outliers else "off", workers=n_workers or 1)
    options.update(station_options)
    result = station(directory_path, lon, lat, **options)
    if debug:
        print(result.status.value_counts())
    accepted = result[result.status.eq("accepted")].copy()
    accepted.attrs["read_errors"] = [issue for entry in result.read_errors[result.read_errors.ne("[]")]
                                      for issue in json.loads(entry)]
    if mnt_corr is not None:
        accepted["wse"] = (accepted.wse - mnt_corr).round(2)
    return accepted


def _process_file(file_path, lon, lat, buffer_size, wse_qual_filter, debug):
    """Single-file adapter for the Observatory incremental cache."""
    from .stations import PRODUCT, extract_candidates
    from pathlib import Path
    path = Path(file_path)
    match = PRODUCT.match(path.name)
    if not match:
        if debug:
            print(f"Unrecognised Raster filename: {path.name}")
        return None
    row = extract_candidates(((path, *match.groups()), []), lon=lon, lat=lat, width_m=None,
                          radius=buffer_size, levels=wse_qual_filter, quality_mode="centre",
                          min_fraction=0, estimator="mean", optical=None, detail_dates=())
    if row["status"] != "accepted":
        return None
    row["wse"] = row["wse_raw"]
    return {key: row[key] for key in ("filename", "date", "cycle", "pass", "resolution", "tile", "wse")}


def extract_qual_timeseries(directory_path, lon, lat, 
                          filter_pass=None, filter_resolution=None, filter_tile=None,  
                          buffer_size=0, debug=False):
    """
    Extract a WSE quality time series at a given position (or around it with buffer)
    from SWOT NetCDF files in a directory (and its subdirectories),
    with options to filter by filename and possibility to remove outliers.

    Parameters:
    -----------
    directory_path : str
        Path to the directory containing SWOT NetCDF files
    lon : float
        Longitude of the point of interest
    lat : float
        Latitude of the point of interest
    filter_pass : str, optional
        Filter by filename (e.g., '394')
    filter_resolution : str, optional
        Filter by filename (e.g., '100')
    buffer_size : int, optional
        Size of the buffer around the point (0 = no buffer, 1 = average of 3x3 pixels, etc.)
    filter_outliers : bool, optional
        If True, remove outliers using IQR method
    debug : bool, optional
        Enable debug output
        
    Returns:
    --------
    pandas.DataFrame
        DataFrame containing dates, qual values, and file paths
    """

    if not os.path.exists(directory_path):
        raise ValueError(f"Directory '{directory_path}' does not exist.")
    
    dates, qual_values, filenames = [], [], []
    passes, resolutions, tiles = [], [], []

    # Find NetCDF files
    nc_files = []
    for root, _, files in os.walk(directory_path):
        for file in files:
            if file.endswith('.nc'):
                if filter_pass and filter_pass not in file:
                    continue
                if filter_resolution and filter_resolution not in file:
                    continue
                if filter_tile and filter_tile not in file:
                    continue
                nc_files.append(os.path.join(root, file))
    
    if not nc_files:
        print(f"No matching NetCDF files found in '{directory_path}'")
        return pd.DataFrame()

    print(f"Processing {len(nc_files)} NetCDF files...")

    for file_path in tqdm(nc_files, desc="Processing files"):
        filename = os.path.basename(file_path)
        
        # Extract pass and resolution from filename
        pass_match = re.search(r'(\d+)', filename)
        res_match = re.search(r'(\d+)m', filename)
        tile_match = re.search(r'(\d+)', filename)
        
        current_pass = pass_match.group(1) if pass_match else "Unknown"
        current_res = res_match.group(1) if res_match else "Unknown"
        current_tile = tile_match.group(1) if tile_match else "Unknown"
        
        if debug:
            print(f"\nProcessing file: {filename}")
            print(f"Pass: {current_pass}, Resolution: {current_res}, Tile: {current_tile}")
        
        try:
            # Extract date from filename
            date_pattern = r'(\d{8}T\d{6})'
            date_matches = re.findall(date_pattern, filename)
            if date_matches:
                date = datetime.strptime(date_matches[0], '%Y%m%dT%H%M%S')
            else:
                date = datetime.fromtimestamp(os.path.getmtime(file_path))
                print(f"No date found in '{filename}', using modification date.")
            
            # Load the file
            swot_file = sload.SWOT_loader(file_path)
            
            if not swot_file.contains(lat, lon):
                if debug:
                    print(f"Skipping file {filename}: point outside bounds")
                swot_file.close()
                continue
            

            # Extract qual value
            try:
                if buffer_size > 0:
                    # Extract mean value in buffer area
                    qual_value = swot_file.extract_area_mean('wse_qual', lon, lat, buffer_size)
                else:
                    # Extract value at single point
                    qual_value = swot_file.extract_point_val('wse_qual', lon, lat)
                
                # Add to results if not NaN
                if not np.isnan(qual_value):
                    dates.append(date)
                    qual_values.append(qual_value)
                    filenames.append(file_path)
                    passes.append(current_pass)
                    resolutions.append(current_res)
                    tiles.append(current_tile)
                elif debug:
                    print(f"Skipping NaN value in {filename}")
                
            except Exception as e:
                print(f"Error extracting WSE from '{filename}': {str(e)}")
                
            swot_file.close()

        except Exception as e:
            print(f"Error processing '{filename}': {str(e)}")

    # Create DataFrame
    df = pd.DataFrame({
        'date': dates,
        'qual': qual_values,
        'filename': filenames,
        'pass': passes,
        'resolution': resolutions,
        'tile': tiles
    })
    
    # Sort by date
    if not df.empty:
        df = df.sort_values(by='date')
            
    return df

def _process_qual_file(file_path, lon, lat, buffer_size, debug):
    """
    Process a single NetCDF file to extract WSE quality data.
    Must be defined at module level for multiprocessing compatibility.
    """
    filename = os.path.basename(file_path)
    
    # Extract pass and resolution from filename
    pass_match = re.search(r'(\d+)', filename)
    res_match = re.search(r'(\d+)m', filename)
    tile_match = re.search(r'(\d+)', filename)
    
    current_pass = pass_match.group(1) if pass_match else "Unknown"
    current_res = res_match.group(1) if res_match else "Unknown"
    current_tile = tile_match.group(1) if tile_match else "Unknown"
    
    if debug:
        print(f"\nProcessing file: {filename}")
        print(f"Pass: {current_pass}, Resolution: {current_res}, Tile: {current_tile}")
    
    try:
        # Extract date from filename
        date_pattern = r'(\d{8}T\d{6})'
        date_matches = re.findall(date_pattern, filename)
        if date_matches:
            date = datetime.strptime(date_matches[0], '%Y%m%dT%H%M%S')
        else:
            date = datetime.fromtimestamp(os.path.getmtime(file_path))
            if debug:
                print(f"No date found in '{filename}', using modification date.")
        
        # Load the file
        swot_file = sload.SWOT_loader(file_path)
        
        if not swot_file.contains(lat, lon):
            if debug:
                print(f"Skipping file {filename}: point outside bounds")
            swot_file.close()
            return None
        
        # Extract qual value
        try:
            if buffer_size > 0:
                # Extract mean value in buffer area
                qual_value = swot_file.extract_area_mean('wse_qual', lon, lat, buffer_size)
            else:
                # Extract value at single point
                qual_value = swot_file.extract_point_val('wse_qual', lon, lat)
            
            swot_file.close()
            
            # Add to results if not NaN
            if not np.isnan(qual_value):
                return {
                    'date': date,
                    'qual': qual_value,
                    'filename': file_path,
                    'pass': current_pass,
                    'resolution': current_res,
                    'tile': current_tile
                }
            elif debug:
                print(f"Skipping NaN value in {filename}")
                return None
            else:
                return None
            
        except Exception as e:
            if debug:
                print(f"Error extracting WSE quality from '{filename}': {str(e)}")
            swot_file.close()
            return None

    except Exception as e:
        if debug:
            print(f"Error processing '{filename}': {str(e)}")
        return None

def extract_qual_timeseries_parallel(directory_path, lon, lat, 
                          filter_pass=None, filter_resolution=None, filter_tile=None,  
                          buffer_size=0, debug=False, n_workers=None):
    """
    Extract a WSE quality time series at a given position (or around it with buffer)
    from SWOT NetCDF files in a directory (and its subdirectories),
    with options to filter by filename and possibility to remove outliers.
    Uses parallel processing to speed up file processing.

    Parameters:
    -----------
    directory_path : str
        Path to the directory containing SWOT NetCDF files
    lon : float
        Longitude of the point of interest
    lat : float
        Latitude of the point of interest
    filter_pass : str, optional
        Filter by filename (e.g., '394')
    filter_resolution : str, optional
        Filter by filename (e.g., '100')
    filter_tile : str, optional
        Filter by tile identifier in filename
    buffer_size : int, optional
        Size of the buffer around the point (0 = no buffer, 1 = average of 3x3 pixels, etc.)
    debug : bool, optional
        Enable debug output
    n_workers : int, optional
        Number of worker processes to use. If None, uses available CPU cores.
        
    Returns:
    --------
    pandas.DataFrame
        DataFrame containing dates, qual values, and file paths
    """
    if not os.path.exists(directory_path):
        raise ValueError(f"Directory '{directory_path}' does not exist.")
    
    # Find NetCDF files
    nc_files = []
    for root, _, files in os.walk(directory_path):
        for file in files:
            if file.endswith('.nc'):
                if filter_pass and filter_pass not in file:
                    continue
                if filter_resolution and filter_resolution not in file:
                    continue
                if filter_tile and filter_tile not in file:
                    continue
                nc_files.append(os.path.join(root, file))
    
    if not nc_files:
        print(f"No matching NetCDF files found in '{directory_path}'")
        # Return empty DataFrame with required columns
        return pd.DataFrame(columns=['date', 'qual', 'filename', 'pass', 'resolution', 'tile'])

    # Set number of workers
    if n_workers is None:
        n_workers = max(1, multiprocessing.cpu_count() - 1)  # Leave one core free
    
    print(f"Processing {len(nc_files)} NetCDF files using {n_workers} workers...")
    
    # Create a partial function with fixed parameters
    process_func = partial(_process_qual_file, 
                          lon=lon, 
                          lat=lat, 
                          buffer_size=buffer_size, 
                          debug=debug)
    
    results = []
    
    # Process files in parallel
    with ProcessPoolExecutor(max_workers=n_workers) as executor:
        # Submit all tasks
        futures = [executor.submit(process_func, file_path) for file_path in nc_files]
        
        # Collect results as they complete
        for future in tqdm(as_completed(futures), total=len(nc_files), desc="Processing files"):
            result = future.result()
            if result is not None:
                results.append(result)

    # Create DataFrame
    if results:
        df = pd.DataFrame(results)
        # Ensure dates are datetime objects (not strings)
        if 'date' in df.columns and not pd.api.types.is_datetime64_any_dtype(df['date']):
            df['date'] = pd.to_datetime(df['date'])
        # Sort by date
        df = df.sort_values(by='date')
    else:
        # Create an empty DataFrame with all required columns
        df = pd.DataFrame(columns=['date', 'qual', 'filename', 'pass', 'resolution', 'tile'])
            
    return df


