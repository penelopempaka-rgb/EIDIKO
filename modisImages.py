import ee, os, requests
import pandas as pd

def runMODIS(polygon_coords, start_date, end_date):
    area = ee.Geometry.Polygon([polygon_coords])
    roi = area.buffer(1000)

    print(f"ROI Area for MODIS: {roi.area().divide(1e6).getInfo():.3f} km²")

    output_dir = os.path.join(os.getcwd(), "data", "MODIS")
    os.makedirs(output_dir, exist_ok=True)

    #MODIS scaling
    def applyModisScaling(image):
        modis_bands = ['ET', 'LE', 'PET', 'ET_QC']
        scaled = image.select(modis_bands).multiply(0.1)
        return image.addBands(scaled, overwrite=True)

    #MODIS collection
    rawModis = (ee.ImageCollection("MODIS/061/MOD16A2")
        .filterBounds(roi)
        .filterDate(start_date, end_date))

    #Filter an image according to the date
    def add_date_property(img):
        return img.set('date_only', img.date().format('YYYY-MM-DD'))

    rawModis = (ee.ImageCollection("MODIS/061/MOD16A2")
        .filterBounds(roi)
        .filterDate(start_date, end_date)
        .map(add_date_property))

    distinct_dates = rawModis.aggregate_array('date_only').distinct()

    def get_first_per_day(d):
        return rawModis.filter(ee.Filter.eq('date_only', d)).first()

    modis = ee.ImageCollection(distinct_dates.map(get_first_per_day)).map(applyModisScaling)

    bands_to_download = ['ET', 'LE', 'PET', 'ET_QC']
    count = modis.size().getInfo()

    if count == 0:
        print("No MODIS images found for this period.")
        return

    print(f"Found {count} MODIS 8-day composites.")
    images_list = modis.toList(count)

    def download_band(img, band, date, out_dir):
        band_img = img.select(band)

        url = band_img.getDownloadURL({
            'scale': 500,
            'region': roi.getInfo(),
            'format': 'GEO_TIFF'
        })

        os.makedirs(out_dir, exist_ok=True)
        out_path = os.path.join(out_dir, f"{band}.tif")

        print(f"Downloading MODIS {date} {band}...")
        r = requests.get(url, stream=True)
        if r.status_code == 200:
            with open(out_path, "wb") as f:
                for chunk in r:
                    f.write(chunk)

    for i in range(count):
        current_img = ee.Image(images_list.get(i))
        timestamp = current_img.get('system:time_start').getInfo()
        date = pd.to_datetime(timestamp, unit='ms').strftime('%Y-%m-%d')

        folder = os.path.join(output_dir, date)
        for band in bands_to_download:
            download_band(current_img, band, date, folder)

    print("Finished downloading all MODIS images.")

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
    
    runMODIS(correctedCoordinates, '2022-01-01', '2023-12-31')
