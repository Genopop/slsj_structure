"""
Download geographic data for SLSJ map generation.

Usage:
    python download_geodata.py
    
This will create a geodata.zip file containing all necessary shapefiles.
"""

import geopandas as gpd
import zipfile
import os
import time
import requests
from pathlib import Path
from shapely.geometry import LineString, Polygon

def download_natural_earth_data():
    """Download required Natural Earth datasets."""
    print("Downloading Natural Earth data...")
    base_url = "https://naturalearth.s3.amazonaws.com"
    
    datasets = {
        'admin': f"{base_url}/50m_cultural/ne_50m_admin_1_states_provinces.zip",
        'lakes': f"{base_url}/10m_physical/ne_10m_lakes.zip",
        'rivers': f"{base_url}/10m_physical/ne_10m_rivers_lake_centerlines.zip",
        'ocean': f"{base_url}/10m_physical/ne_10m_ocean.zip"
    }
    
    gdfs = {}
    for name, url in datasets.items():
        print(f"  Downloading {name}...")
        gdfs[name] = gpd.read_file(url)
        print(f"    ✓ {name}: {len(gdfs[name])} features")
    
    return gdfs

def download_osm_rivers():
    """Download river data from OpenStreetMap for SLSJ region."""
    print("\nDownloading OpenStreetMap rivers for SLSJ region...")
    
    # Bounding box for SLSJ region (with buffer)
    min_lat, min_lon = 47.9, -72.8
    max_lat, max_lon = 49.2, -70.0
    
    buffer = 0.8
    bbox = (min_lat - buffer, min_lon - buffer, max_lat + buffer, max_lon + buffer)
    
    overpass_url = "http://overpass-api.de/api/interpreter"
    
    # Match the original script's filtered query for specific rivers only
    query = f"""
    [out:json][timeout:180];
    (
      way["waterway"="river"]["name"~"Péribonka|Métabetchouane|Belle Rivière|Mistassini|Ashuapmushuan|Ha! Ha!|Saguenay|Décharge",i]({bbox[0]},{bbox[1]},{bbox[2]},{bbox[3]});
      relation["waterway"="river"]["name"~"Péribonka|Métabetchouane|Belle Rivière|Mistassini|Ashuapmushuan|Ha! Ha!|Saguenay|Décharge",i]({bbox[0]},{bbox[1]},{bbox[2]},{bbox[3]});
      way["waterway"="riverbank"]["name"~"Saguenay|Ha! Ha!|Décharge",i]({bbox[0]},{bbox[1]},{bbox[2]},{bbox[3]});
      relation["waterway"="riverbank"]["name"~"Saguenay|Ha! Ha!|Décharge",i]({bbox[0]},{bbox[1]},{bbox[2]},{bbox[3]});
      way["natural"="water"]["name"~"Saguenay|Ha! Ha!|Décharge",i]({bbox[0]},{bbox[1]},{bbox[2]},{bbox[3]});
      relation["natural"="water"]["name"~"Saguenay|Ha! Ha!|Décharge",i]({bbox[0]},{bbox[1]},{bbox[2]},{bbox[3]});
      way["natural"="bay"]["name"~"Ha.?Ha",i]({bbox[0]},{bbox[1]},{bbox[2]},{bbox[3]});
    );
    (._;>;);
    out geom;
    """
    
    max_retries = 3
    for attempt in range(max_retries):
        try:
            print(f"  Querying Overpass API (attempt {attempt + 1}/{max_retries})...")
            response = requests.post(overpass_url, data={'data': query}, timeout=180)
            response.raise_for_status()
            data = response.json()
            break
        except Exception as e:
            if attempt < max_retries - 1:
                print(f"    Failed: {e}. Retrying in 2 seconds...")
                time.sleep(2)
            else:
                print(f"    All attempts failed. OSM data will be missing.")
                return None, None
    
    # Process data
    river_lines = []
    river_names = []
    water_polygons = []
    
    for element in data.get('elements', []):
        if element['type'] == 'way' and 'geometry' in element:
            coords = [(node['lon'], node['lat']) for node in element.get('geometry', [])]
            if len(coords) < 2:
                continue
            
            tags = element.get('tags', {})
            name = tags.get('name', 'Unknown')
            
            # Check if it's an area (riverbank, water, bay)
            is_area = (tags.get('waterway') == 'riverbank' or 
                      tags.get('natural') in ['water', 'bay'])
            
            if is_area:
                if len(coords) >= 3 and coords[0] == coords[-1]:
                    water_polygons.append(Polygon(coords))
            else:
                # It's a river centerline
                river_lines.append(LineString(coords))
                river_names.append(name)
    
    print(f"    ✓ Retrieved {len(river_lines)} river segments")
    print(f"    ✓ Retrieved {len(water_polygons)} water polygons")
    
    # Create GeoDataFrames
    osm_rivers = None
    osm_water = None
    
    if river_lines:
        osm_rivers = gpd.GeoDataFrame(
            {'name': river_names, 'geometry': river_lines},
            crs='EPSG:4326'
        )
    
    if water_polygons:
        osm_water = gpd.GeoDataFrame(
            {'geometry': water_polygons},
            crs='EPSG:4326'
        )
    
    return osm_rivers, osm_water

def save_to_zip(gdfs, output_path='geodata.zip'):
    """Save GeoDataFrames to a zip file."""
    print(f"\nSaving to {output_path}...")
    
    # Create temporary directory
    temp_dir = Path('temp_geodata')
    temp_dir.mkdir(exist_ok=True)
    
    try:
        # Save each GeoDataFrame as shapefile
        for name, gdf in gdfs.items():
            print(f"  Saving {name}...")
            # Create subdirectory for this layer
            layer_dir = temp_dir / name
            layer_dir.mkdir(exist_ok=True)
            gdf.to_file(layer_dir / f'{name}.shp')
        
        # Create zip file
        with zipfile.ZipFile(output_path, 'w', zipfile.ZIP_DEFLATED) as zipf:
            for name in gdfs.keys():
                folder = temp_dir / name
                for file in folder.rglob('*'):
                    if file.is_file():
                        arcname = f"{name}/{file.name}"
                        zipf.write(file, arcname)
                        
        print(f"✓ Created {output_path}")
        
        # Get file size
        size_mb = os.path.getsize(output_path) / (1024 * 1024)
        print(f"  File size: {size_mb:.1f} MB")
        
    finally:
        # Clean up temp directory
        import shutil
        if temp_dir.exists():
            shutil.rmtree(temp_dir)

def main():
    print("=" * 60)
    print("Geographic Data Download Script for SLSJ Map")
    print("=" * 60)
    print()
    
    # Download Natural Earth data
    gdfs = download_natural_earth_data()
    
    # Download OSM data
    osm_rivers, osm_water = download_osm_rivers()
    if osm_rivers is not None:
        gdfs['osm_rivers'] = osm_rivers
    if osm_water is not None:
        gdfs['osm_water'] = osm_water
    
    # Save to zip
    save_to_zip(gdfs)
    
    print()
    print("=" * 60)
    print("SUCCESS!")
    print("=" * 60)
    print()

if __name__ == "__main__":
    main()