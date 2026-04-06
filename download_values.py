import os
import glob
import rasterio
import pandas as pd
import numpy as np
from scipy import stats
from datetime import datetime, timedelta

BASE_DIR = os.getcwd()
LANDSAT_DIR = os.path.join(BASE_DIR, "data", "LANDSAT")
MODIS_DIR = os.path.join(BASE_DIR, "data", "MODIS")
OUTPUT_CSV = os.path.join(BASE_DIR, "et_training_data.csv")


def calculate_dry_edges(ndvi, lst):
    #Filet values. We ignore NaN and non logical NDVI values
    mask = (~np.isnan(ndvi)) & (~np.isnan(lst)) & (ndvi >= -1) & (ndvi <= 1)
    ndvi_values = ndvi[mask]
    lst_values = lst[mask]

    #Seperate NDVI in "bins"
    bins = np.arange(0, 1.05, 0.05)
    max_list_per_bin = []
    ndvi_centers = []

    for i in range(len(bins) - 1):
        #find pixels that belong to the current bin
        index = (ndvi_values >= bins[i]) & (ndvi_values < bins[i + 1])

        if np.any(index):
            #keep 99% of temperature to avoid outliers
            max_t = np.percentile(lst_values[index], 99)
            max_list_per_bin.append(max_t)
            ndvi_centers.append((bins[i] + bins[i+1]) / 2)

    #Perform linear regression to find a(intercept) and b(slope)
    if len(max_list_per_bin) > 2:
        b, a, r_value, p_value, std_err = stats.linregress(ndvi_centers, max_list_per_bin)
        return a, b
    else:
        return None, None

def calculate(bands):
    c = 0.5569
    x = 1e-10           #to avoid division by 0
    results = {}

    r1 = bands['B1'] / 10000   # Ultra blue
    r2 = bands['B2'] / 10000   # Blue  
    r3 = bands['B3'] / 10000   # Green   
    r4 = bands['B4'] / 10000   # Red   
    r5 = bands['B5'] / 10000   # Surface Reflectance (near infrared)   
    r6 = bands['B6'] / 10000   # Surface Reflectance (shortwave infrared 1)   
    r7 = bands['B7'] / 10000   # Surface Reflectance (shortwave infrared 2)   

    #Greeness Indices
    results['NDVI'] = (r5 - r4)/(r5 + r4 + x)
    results['EVI'] = 2.5 * (r5 - r4) / (r5 + 6 * r4 - 7.5 * r2 + 1 + x)
    results['SAVI']= 1.5 * ((r5 - r4) / (r5 + r4 + 0.5 + x))

    #MSAVI calculation
    results['MSAVI'] = (2 * r5 + 1 - np.sqrt((2 * r5 + 1)**2 - 8*(r5 - r4))) / 2

    #Vegetation Water Indices
    results['NDMI']= (r5 - r6) / (r5 + r6) + x
    results['NDWI'] = (r3 - r5) / (r3 + r5 + x)
    results['NDIIr7'] = (r5 - r7) / (r5 + r7 + x)
    results['D1609'] = 1 - (r6 / ((1 - c) * r5 + c * r7 + x))

    #Aldebo using weigths from Table 4
    results['Albedo'] = r1 * 0.130 + r2 * 0.115 + r3 * 0.143 + r4 * 0.281 + r5 * 0.180 + r6 * 0.108 + r7 * 0.042

    #Land Surface Temperature (LST)
    results['LST'] = bands['B10'] 

    #Temperature Vegetation Dryness Index (TVDI)
    ts_min = np.nanmin(results['LST'])
    a, b = calculate_dry_edges(results['NDVI'], results['LST'])
    if a is not None and b is not None:
        dry_edges = a + b * results['NDVI']
        results['TDVI'] = (results['LST'] - ts_min) / (dry_edges - ts_min + x)
    else:
        results['TDVI'] = (results['LST'] - ts_min) / (np.nanmax(results['LST']) - ts_min)

    return results

def main():
    landsat_folders = sorted(glob.glob(os.path.join(LANDSAT_DIR, "*")))

    if not landsat_folders:
        print("No folders in LANDSAT DIR.")
        return
    
    for folder in landsat_folders:
        date_str = os.path.basename(folder)
        print(f"Date proccessing: {date_str}")
        try:
            bands_paths = {
                'B1': glob.glob(os.path.join(folder, "*B1.tif"))[0],
                'B2': glob.glob(os.path.join(folder, "*B2.tif"))[0],
                'B3': glob.glob(os.path.join(folder, "*B3.tif"))[0],
                'B4': glob.glob(os.path.join(folder, "*B4.tif"))[0],
                'B5': glob.glob(os.path.join(folder, "*B5.tif"))[0],
                'B6': glob.glob(os.path.join(folder, "*B6.tif"))[0],
                'B7': glob.glob(os.path.join(folder, "*B7.tif"))[0],
                'B10': glob.glob(os.path.join(folder, "*B10.tif"))[0],
            }
        except IndexError:
            print(f"Error: Missing bands in folder {date_str}")
            continue

        bands_data = {}

        with rasterio.open(bands_paths['B4']) as ref_src:
            target_shape = (ref_src.height, ref_src.width)

        for b_name, b_path in bands_paths.items():
            with rasterio.open(b_path) as src:
                #Force all bands (like B10) to match the B4 shape using bilinear resampling
                bands_data[b_name] = src.read(
                    1,
                    out_shape = target_shape,
                    resampling = rasterio.enums.Resampling.bilinear
                ).astype('float32')

        results = calculate(bands_data)

        #Upscaling to 1km (Aggregation)
        upscaled_data = {}
        window_size = 33
        total_pixels_per_grid = window_size * window_size

        for name, arr in results.items():
            #Adjust dimensions so that they are myltiples of 33
            new_h = arr.shape[0] // window_size
            new_w = arr.shape[1] // window_size
            temp = arr[:new_h * window_size, :new_w * window_size]


            #Reorganization in blocks
            reshaped = temp.reshape(new_h, window_size, new_w, window_size)

            #Calculate NaN percentage. Count how many NaN pixels per block
            nan_count = np.isnan(reshaped).sum(axis=(1, 3))
            nan_percentage = nan_count / total_pixels_per_grid

            #Calculate mea value (Aggregation)
            grid_mean = np.nanmean(reshaped, axis=(1, 3))

            #If clouds > 30% set the whole 1km grid as NaN
            grid_mean[nan_percentage > 0.50] = np.nan

            upscaled_data[name] = grid_mean

        landsat_dt = datetime.strptime(date_str, '%Y-%m-%d')
        modis_path = None
        matched_modis_date = None

        #Search MODIS that starts 1-4 days AGO
        for offset in range(0, 9):
            check_date = (landsat_dt - timedelta(days=offset)).strftime('%Y-%m-%d')
            potential_path = os.path.join(MODIS_DIR, check_date, "ET.tif")

            if os.path.exists(potential_path):
                modis_path = potential_path
                matched_modis_date = check_date
                break

        #If no pictures found, skip
        if not modis_path:
            print(f"Skipping {date_str}: No MODIS ET found in 1-4 day window")
            continue

            print(f"Match found: LANDSAT {date_str} -> MODIS {matched_modis_date}")

        if not os.path.exists(modis_path):
            print(f"No MODIS ET in the folder {date_str}")
            continue

        with rasterio.open(modis_path) as modis_source:
            target_h, target_w = next(iter(upscaled_data.values())).shape

            #Resample MODIS ET to match the Landsat aggregated grid
            modis_et = modis_source.read(
                1,
                out_shape = (target_h, target_w),
                resampling = rasterio.enums.Resampling.bilinear
            ).astype('float32')

            modis_et = modis_et * 0.1

        df_dict = {name: arr.flatten() for name, arr in upscaled_data.items()}
        df_dict['MODIS_ET'] = modis_et.flatten()

        #Convert each image to data (DataFrame and Flattening)
        df_day = pd.DataFrame(df_dict)

        #Remove NaN values from clouds
        df_day = df_day.dropna()
        df_day['date'] = date_str

        #Save on CSV
        if not os.path.exists(OUTPUT_CSV):
            df_day.to_csv(OUTPUT_CSV, index=False)
        else:
            df_day.to_csv(OUTPUT_CSV, mode = 'a', header=False, index=False)

        print(f"Date {date_str} is complete. {len(df_day)} were added.")

if __name__ == "__main__":
    main()

