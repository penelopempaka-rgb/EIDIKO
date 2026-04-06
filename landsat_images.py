import ee, os, requests
import argparse
import pandas as pd
import shutil

def runLandsat(polygon_coords, start_date, end_date):
    
    area = ee.Geometry.Polygon([polygon_coords])
    
    print("Area size:", area.area().divide(1e6).getInfo(), "km²")
    roi = area.buffer(1000)
    roi_area_m2 = roi.area().getInfo()
    roi_area_km2 = roi_area_m2 / 1e6
    print(f"ROI Area: {roi_area_km2:.3f} km²")
    
    output_dir = os.path.join(os.getcwd(), "data", "LANDSAT")
    # Clear folder before each Landsat 
    # if os.path.exists(output_dir):
    #     shutil.rmtree(output_dir)
    # os.makedirs(output_dir, exist_ok=True)
    
    def cloudMask(image):
      # Define cloud shadow and cloud bitmasks (Bits 3 and 5)
      cloudShadowBitmask = (1 << 3)
      cloudBitmask = (1 << 5)
    
      # Select the Quality Assessment (QA) band for pixel quality information
      qa = image.select('QA_PIXEL')
    
      # Create a binary mask to identify clear conditions (both cloud and cloud shadow bits set to 0)
      mask = (qa.bitwiseAnd(cloudShadowBitmask).eq(0)
                .And(qa.bitwiseAnd(cloudBitmask).eq(0)))
    
      # Update the original image, masking out cloud and cloud shadow-affected pixels
      return image.updateMask(mask)
    
    
    def applyScaleFactors(image):
        # Surface Reflectance bands
        sr_bands = ['SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B6', 'SR_B7']
        # Thermal and related bands
        thermal_band = 'ST_B10'
        emissivity_band = 'ST_EMIS'
        atrans_band = 'ST_ATRAN'
    
        # Scale surface reflectance bands: multiply by 0.0000275 and subtract 0.2
        scaled_sr = image.select(sr_bands).multiply(0.0000275).add(-0.2)
    
        # Scale thermal band: multiply by 0.00341802 and add 149.0
        scaled_thermal = image.select(thermal_band).multiply(0.00341802).add(149.0)
    
        # Scale emissivity and atmospheric transmissivity bands
        scaled_emis = image.select(emissivity_band).multiply(0.0001)
        scaled_atran = image.select(atrans_band).multiply(0.0001)
    
        # Replace old bands with scaled bands
        image = image.addBands(scaled_sr, overwrite=True)
        image = image.addBands(scaled_thermal, overwrite=True)
        image = image.addBands(scaled_emis, overwrite=True)
        image = image.addBands(scaled_atran, overwrite=True)
        return image

    # --- START: Print ust to know all the raw data
    # landsat_raw = ee.ImageCollection("LANDSAT/LC08/C02/T1_L2") \
    # .filterBounds(roi) \
    # .filterDate(start_date, end_date) \
    # .sort("system:time_start", False)

    # collection_raw = landsat_raw.toList(landsat_raw.size())
    # # Loop through them
    # count = collection_raw.size().getInfo()
    # print(f"Found {count} Landsat images")
    # for i in range(count):
    #     img = collection_raw.get(i)
    #     img = ee.Image(img)
    #     time_start = img.get('system:time_start').getInfo()
    #     date_str = datetime.datetime.utcfromtimestamp(time_start / 1000).strftime('%Y-%m-%d')
    #     print(f"Image {i+1}: {date_str}")
    # ---- END
    
#     # Use Landsat Satellite to create the collection and keep the number of images found
#     landsat = (ee.ImageCollection("LANDSAT/LC08/C02/T1_L2")
#         .filterBounds(roi)
#         .filterDate(start_date, end_date)
#         .filter(ee.Filter.lt("CLOUD_COVER", 60))
#         .map(applyScaleFactors))
# #        .map(cloudMask))   # Removed cloud masking because it removes pixels when cloudiness ~ 60% and fails to calculate hot/cold pixels
    
    rawLandsat = (ee.ImageCollection("LANDSAT/LC08/C02/T1_L2")
        .filterBounds(roi)
        .filterDate(start_date, end_date)
        .filter(ee.Filter.lt("CLOUD_COVER", 60)))

    # Filtering so that we keep only one image per day
    # Ensure that there are no duplicates
    def add_date_property(img):
        return img.set('date_only', img.date().format('YYYY-MM-DD'))

    rawLandsatWithDate = rawLandsat.map(add_date_property)

    def keep_distinct_days(collection):
        distinct_dates = collection.aggregate_array('date_only').distinct()
        
        def get_first_per_day(d):
            return collection.filter(ee.Filter.eq('date_only', d)).first()
            
        return ee.ImageCollection(distinct_dates.map(get_first_image_per_day))

    def get_first_image_per_day(d):
        return rawLandsatWithDate.filter(ee.Filter.eq('date_only', d)).first()

    landsat = ee.ImageCollection(rawLandsatWithDate.aggregate_array('date_only').distinct().map(get_first_image_per_day)).map(applyScaleFactors)

    bands = ['SR_B1', 'SR_B2','SR_B3','SR_B4','SR_B5','SR_B6','SR_B7','ST_B10','ST_EMIS', 'ST_ATRAN']
    collection_size = landsat.size().getInfo()
    
    if collection_size==0: 
        print("No LANDSAT Images to download the last 15 days")
        return
    else:
        images = landsat.toList(landsat.size())
        count = images.size().getInfo()
        print(f"Found {count} Landsat images.")
    
        # Function to download each band and save into output_dir
        def download_band(img, band, date, out_dir):
            # if band == 'ST_B10':
            #     band_img = img.select(band).updateMask(img.select(band).gt(150))
            # else:
            #     band_img = img.select(band).updateMask(img.select(band).gt(0))
            band_img = img.select(band)
            if band != 'ST_B10':
                band_img = band_img.updateMask(band_img.gte(-0.2))       # keep values above minimal scaled
                
            # Set proper scale in meters
            if band in ['ST_B10', 'ST_EMIS', 'ST_ATRAN']:
                scale = 100   # This is the roght one for Landsat 8
                #scale = 30   # because my roi is super small
            else:
                scale = 30
        
            url = band_img.getDownloadURL({
                'scale': scale,
                'region': roi.getInfo(),
                'format': 'GEO_TIFF'
            })
            
            os.makedirs(out_dir, exist_ok=True)
            out_path = f"{out_dir}/{band}.tif"
        
            print(f"Downloading {date} {band} -> {out_path}")
            r = requests.get(url, stream=True)
            with open(out_path, "wb") as f:
                for chunk in r:
                    f.write(chunk)
        
        # For each image seperate bands using download_band func
        for i in range(count):
            mg = ee.Image(images.get(i))
            #latest_img = landsat.sort("system:time_start", False).first()
            timestamp = mg.get('system:time_start').getInfo()
            date = pd.to_datetime(timestamp, unit='ms').strftime('%Y-%m-%d')
            print(f"Image {i+1} from {count}: {date}")
            folder = f"{output_dir}/{date}"
            os.makedirs(folder, exist_ok=True)
        
            for band in bands:
                download_band(mg, band, date, folder)
        
        # Download DEM once
        def download_dem():
            #dem = ee.Image("USGS/SRTMGL1_003").clip(roi)
            dem = ee.ImageCollection("COPERNICUS/DEM/GLO30").select("DEM").mosaic().clip(roi)
            url = dem.getDownloadURL({
                'scale': 30,
                'region': roi.getInfo(),
                'format': 'GEO_TIFF',
                'maxPixels': 1e13
            })
    
            dem_dir = os.path.join(os.getcwd(), "data", "automate", "DEM")
            os.makedirs(dem_dir, exist_ok=True)
            dem_file = os.path.join(dem_dir, "dem.tif")
            print("Started DEM export task")
            r = requests.get(url, stream=True)
            with open(dem_file, "wb") as f:
                for chunk in r:
                    f.write(chunk)
        
        download_dem()
        
        print("Finished downloading all satellite images")
        return images

if __name__ == "__main__":
    ee.Authenticate()
    ee.Initialize(project = 'ee-penelopempaka')

        #Define the coordinates of the field
    geojson = {
            "type": "Polygon",
            "coordinates": [
                [
                    [39.68549395514354, 22.77384025414311],
                    [39.68542010286271, 22.77388105364505],
                    [39.6853661717447, 22.77395131437986],
                    [39.68532643518604, 22.7740361890915],
                    [39.68530662710225, 22.77411013994296],
                    [39.68529252391261, 22.77416194664032],
                    [39.6852841540361, 22.77423219436491],
                    [39.68530393883208, 22.77435364786365],
                    [39.6855817198441, 22.77516765090019],
                    [39.68621575750824, 22.77478349947132],
                    [39.68578785932531, 22.77367987358766],
                    [39.68549395514354, 22.77384025414311]
                ]
            ]
    }

    #[Lat, Lon] -> [Lon, Lat] for GEE
    rawPoints = geojson['coordinates'][0]
    correctedCoordinates = [[point[1], point[0]] for point in rawPoints]
    
    runLandsat(correctedCoordinates, '2022-01-01', '2023-12-31')