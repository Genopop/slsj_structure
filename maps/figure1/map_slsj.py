"""
Saguenay–Lac-Saint-Jean Population Distribution Map Generator

This script generates a cartographic visualization of population
distribution across the Saguenay–Lac-Saint-Jean region of Quebec, Canada. The map
displays municipalities as proportionally-sized circles colored by watercourse 
subdivision, with geographic features from Natural Earth and OpenStreetMap.

Author: Gilles-Philippe Morin, with the help of Gemini 3 Pro and Claude Sonnet 4.5
"""

import warnings
import time
import math
from typing import Dict, List, Tuple, Optional, Union

import numpy as np
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.patheffects as pe
from shapely.geometry import box, LineString, MultiLineString, Polygon
from shapely.ops import linemerge, unary_union
from scipy.interpolate import make_interp_spline
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import requests

warnings.filterwarnings('ignore')


# =============================================================================
# CONFIGURATION - COLOR SCHEMES
# =============================================================================

SUBDIVISION_COLORS: Dict[str, str] = {
    "East of Ha! Ha!": 'gold',
    "South of Saguenay River": 'lime',
    "North of Saguenay River": 'deeppink',
    "East of Lac Saint-Jean": 'darkgreen',
    "La-Belle-Rivière–Métabetchouane": 'aqua',
    "Mistassini-Péribonka": 'rosybrown',
    "Métabetchouane-Ashuapmushuan": 'dodgerblue',
    "Ashuapmushuan-Mistassini": '#a020f0'
}

CITY_SUBDIVISIONS: Dict[str, str] = {
    "Albanel": "Ashuapmushuan-Mistassini",
    "Alma": "East of Lac Saint-Jean",
    "Arvida": "South of Saguenay River",
    "Bégin": "North of Saguenay River",
    "Chambord": "Métabetchouane-Ashuapmushuan",
    "Chicoutimi": "South of Saguenay River",
    "Chicoutimi-Nord": "North of Saguenay River",
    "Chute-des-Passes": "Mistassini-Péribonka",
    "Delisle": "East of Lac Saint-Jean",
    "Desbiens": "La-Belle-Rivière–Métabetchouane",
    "Dolbeau": "Ashuapmushuan-Mistassini",
    "Ferland-et-Boilleau": "South of Saguenay River",
    "Girardville": "Ashuapmushuan-Mistassini",
    "Hébertville": "East of Lac Saint-Jean",
    "Hébertville-Station": "East of Lac Saint-Jean",
    "Jonquière": "South of Saguenay River",
    "L'Anse-Saint-Jean": "East of Ha! Ha!",
    "L'Ascension-de-Notre-Seigneur": "East of Lac Saint-Jean",
    "La Baie": "South of Saguenay River",
    "La Doré": "Métabetchouane-Ashuapmushuan",
    "Labrecque": "East of Lac Saint-Jean",
    "Lac-Bouchette": "Métabetchouane-Ashuapmushuan",
    "Lac-Kénogami": "South of Saguenay River",
    "Lac-à-la-Croix": "La-Belle-Rivière–Métabetchouane",
    "Larouche": "South of Saguenay River",
    "Laterrière": "South of Saguenay River",
    "Mashteuiatsh": "Métabetchouane-Ashuapmushuan",
    "Mistassini": "Mistassini-Péribonka",
    "Métabetchouan": "La-Belle-Rivière–Métabetchouane",
    "Mont-Apica": "East of Lac Saint-Jean",
    "Normandin": "Ashuapmushuan-Mistassini",
    "Notre-Dame-de-Lorette": "Mistassini-Péribonka",
    "Notre-Dame-du-Rosaire": "East of Lac Saint-Jean",
    "Petit-Saguenay": "East of Ha! Ha!",
    "Péribonka": "Mistassini-Péribonka",
    "Rivière-Éternité": "East of Ha! Ha!",
    "Roberval": "Métabetchouane-Ashuapmushuan",
    "Saint-Ambroise": "North of Saguenay River",
    "Saint-André-du-Lac-Saint-Jean": "La-Belle-Rivière–Métabetchouane",
    "Saint-Augustin-du-Lac-Saint-Jean": "Mistassini-Péribonka",
    "Saint-Bruno": "East of Lac Saint-Jean",
    "Saint-Charles-de-Bourget": "North of Saguenay River",
    "Saint-David-de-Falardeau": "North of Saguenay River",
    "Saint-Edmond": "Ashuapmushuan-Mistassini",
    "Saint-Eugène-du-Lac-Saint-Jean": "Mistassini-Péribonka",
    "Saint-François-de-Sales": "Métabetchouane-Ashuapmushuan",
    "Saint-Fulgence": "North of Saguenay River",
    "Saint-Félicien": "Métabetchouane-Ashuapmushuan",
    "Saint-Félix-d'Otis": "East of Ha! Ha!",
    "Saint-Gédéon": "East of Lac Saint-Jean",
    "Saint-Henri-de-Taillon": "East of Lac Saint-Jean",
    "Saint-Honoré-de-Chicoutimi": "North of Saguenay River",
    "Saint-Ludger-de-Milot": "Mistassini-Péribonka",
    "Saint-Méthode": "Ashuapmushuan-Mistassini",
    "Saint-Nazaire": "East of Lac Saint-Jean",
    "Saint-Prime": "Métabetchouane-Ashuapmushuan",
    "Saint-Stanislas-du-Lac-Saint-Jean": "Mistassini-Péribonka",
    "Saint-Thomas-Didyme": "Ashuapmushuan-Mistassini",
    "Sainte-Hedwidge": "Métabetchouane-Ashuapmushuan",
    "Sainte-Jeanne-d'Arc-du-Lac-Saint-Jean": "Mistassini-Péribonka",
    "Sainte-Monique-de-Honfleur": "East of Lac Saint-Jean",
    "Sainte-Rose-du-Nord": "North of Saguenay River",
    "Sainte-Élisabeth-de-Proulx": "Mistassini-Péribonka",
    "Shipshaw": "North of Saguenay River",
    "Val-Jalbert": "Métabetchouane-Ashuapmushuan"
}


# =============================================================================
# CONFIGURATION - VISUALIZATION PARAMETERS
# =============================================================================

# Coordinate reference systems
TARGET_CRS = "EPSG:32198"  # NAD83 / Quebec Lambert
SOURCE_CRS = "EPSG:4326"   # WGS84

# Colors
WATER_COLOR = '#a3ccff'
WATER_LABEL_COLOR = '#3a6ea5'
RIVER_LABEL_COLOR = '#555555'
BACKGROUND_COLOR = 'white'

# Spatial parameters (in degrees for source data, meters for projected data)
BUFFER_DEGREES = 0.8
BUFFER_METERS = 25000

# Population visualization
POPULATION_SCALE_FACTOR = 40.0
DEFAULT_POPULATION = 1000
FALLBACK_POPULATION = 500
LABEL_SIZE_THRESHOLD = 2500

# Circle dispersion algorithm
DISPERSION_ITERATIONS = 50
REPULSION_FACTOR = 0.5
ATTRACTION_FACTOR = 0.05
PADDING_FACTOR = 1.05

# Text labeling
RIVER_LABEL_MIN_LENGTH = 2000  # Minimum river length (m) for labeling
ANGLE_CALC_DELTA = 100         # Distance (m) for tangent angle calculation
TEXT_PATH_COVERAGE = 0.8       # Fraction of path used for text placement

# Label positioning offsets (in meters)
LAC_SJ_LABEL_OFFSET_X = 3000
LAC_SJ_LABEL_OFFSET_Y_TOP = 5000
LAC_SJ_LABEL_OFFSET_Y_BOTTOM = -3000
SAGUENAY_LABEL_OFFSET_Y_TOP = 5000
SAGUENAY_LABEL_OFFSET_Y_BOTTOM = -3000
SAGUENAY_DEFAULT_OFFSET_Y = 6000
SAGUENAY_RAY_LENGTH = 20000
CITY_LABEL_OFFSET = 500

# Map dimensions and layout
FIGURE_SIZE = (10, 10)
INSET_SIZE_FRACTION = 0.30
INSET_PADDING_FRACTION = 0.1

# Legend parameters
LEGEND_GAP = 0.015
LEGEND_PAD = 0.015
LEGEND_TITLE_HEIGHT = 0.025
LEGEND_LABEL_HEIGHT = 0.020
LEGEND_SPACING = 0.010
LEGEND_POPULATION_SIZES = [100, 1000, 10000]
LEGEND_POPULATION_LABELS = ['100', '1,000', '10,000']

# API parameters
OSM_TIMEOUT = 180
OSM_MAX_RETRIES = 3
OSM_RETRY_DELAY = 2


# =============================================================================
# DATA LOADING
# =============================================================================

def load_population_data(filepath: str) -> gpd.GeoDataFrame:
    """
    Load population data and convert to GeoDataFrame.
    
    Parameters
    ----------
    filepath : str
        Path to CSV file containing population data with columns:
        'Longitude', 'Latitude', "Location of proband's marriage", 'Population size'
    
    Returns
    -------
    gpd.GeoDataFrame
        GeoDataFrame with point geometries and mapped subdivision colors
    """
    df = pd.read_csv(filepath)
    df['Color'] = (df["Location of proband's marriage"]
                   .map(CITY_SUBDIVISIONS)
                   .map(SUBDIVISION_COLORS))
    df['Population size'] = df.get('Population size', DEFAULT_POPULATION).fillna(FALLBACK_POPULATION)
    
    return gpd.GeoDataFrame(
        df,
        geometry=gpd.points_from_xy(df.Longitude, df.Latitude),
        crs=SOURCE_CRS
    )


def download_natural_earth_data() -> Tuple[gpd.GeoDataFrame, gpd.GeoDataFrame, 
                                           gpd.GeoDataFrame, gpd.GeoDataFrame]:
    """
    Download required Natural Earth datasets.
    
    Returns
    -------
    tuple of gpd.GeoDataFrame
        (admin_boundaries, lakes, rivers, ocean) datasets
    """
    print("Downloading Natural Earth data...")
    base_url = "https://naturalearth.s3.amazonaws.com"
    urls = [
        f"{base_url}/50m_cultural/ne_50m_admin_1_states_provinces.zip",
        f"{base_url}/10m_physical/ne_10m_lakes.zip",
        f"{base_url}/10m_physical/ne_10m_rivers_lake_centerlines.zip",
        f"{base_url}/10m_physical/ne_10m_ocean.zip"
    ]
    return tuple(gpd.read_file(url) for url in urls)


def query_openstreetmap_rivers(bbox: Tuple[float, float, float, float]) -> Tuple[Dict[str, List[LineString]], List[Polygon]]:
    """
    Query OpenStreetMap Overpass API for river data.
    
    Parameters
    ----------
    bbox : tuple of float
        Bounding box as (min_lat, min_lon, max_lat, max_lon) in WGS84
    
    Returns
    -------
    tuple
        (osm_lines, osm_polygons) where:
        - osm_lines: dict mapping river names to lists of LineString geometries
        - osm_polygons: list of Polygon geometries for water bodies
    """
    print("Querying OpenStreetMap...")
    
    query = f"""
    [out:json][timeout:{OSM_TIMEOUT}];
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
    
    osm_lines: Dict[str, List[LineString]] = {}
    osm_polygons: List[Polygon] = []
    
    for attempt in range(OSM_MAX_RETRIES):
        try:
            response = requests.get(
                "https://overpass-api.de/api/interpreter",
                params={'data': query},
                headers={'User-Agent': 'PythonMapScript/18.0'},
                timeout=OSM_TIMEOUT + 10
            )
            
            if response.status_code == 200:
                for elem in response.json().get('elements', []):
                    if elem['type'] != 'way' or 'geometry' not in elem:
                        continue
                    
                    tags = elem.get('tags', {})
                    name = tags.get('name', 'Unknown')
                    coords = [(p['lon'], p['lat']) for p in elem['geometry']]
                    
                    if len(coords) < 2:
                        continue
                    
                    is_area = (tags.get('waterway') == 'riverbank' or 
                             tags.get('natural') in ['water', 'bay'])
                    
                    if is_area:
                        osm_polygons.append(Polygon(coords))
                    else:
                        osm_lines.setdefault(name, []).append(LineString(coords))
                break
        except (requests.RequestException, ValueError) as e:
            if attempt < OSM_MAX_RETRIES - 1:
                time.sleep(OSM_RETRY_DELAY)
            else:
                print(f"Warning: OSM query failed after {OSM_MAX_RETRIES} attempts: {e}")
    
    return osm_lines, osm_polygons


# =============================================================================
# GEOMETRY PROCESSING
# =============================================================================

def project_and_clip_layers(
    gdf: gpd.GeoDataFrame,
    lakes_gdf: gpd.GeoDataFrame,
    rivers_gdf: gpd.GeoDataFrame,
    ocean_gdf: gpd.GeoDataFrame,
    osm_lines: Dict[str, List[LineString]],
    osm_polygons: List[Polygon]
) -> Tuple[gpd.GeoDataFrame, Polygon, gpd.GeoDataFrame, gpd.GeoDataFrame,
           gpd.GeoDataFrame, gpd.GeoDataFrame, gpd.GeoDataFrame, List[Polygon]]:
    """
    Project all layers to target CRS and clip to viewport.
    
    Parameters
    ----------
    gdf : gpd.GeoDataFrame
        Population point data in WGS84
    lakes_gdf, rivers_gdf, ocean_gdf : gpd.GeoDataFrame
        Natural Earth layers in WGS84
    osm_lines : dict
        OSM river centerlines by name
    osm_polygons : list
        OSM water body polygons
    
    Returns
    -------
    tuple
        (gdf_proj, viewport_box, viewport_gdf, lakes_final, rivers_final,
         rivers_proj, ocean_final, water_polygons)
    """
    gdf_proj = gdf.to_crs(TARGET_CRS)
    minx, miny, maxx, maxy = gdf_proj.total_bounds
    
    viewport_box = box(
        minx - BUFFER_METERS, miny - BUFFER_METERS,
        maxx + BUFFER_METERS, maxy + BUFFER_METERS
    )
    viewport_gdf = gpd.GeoDataFrame({'geometry': [viewport_box]}, crs=TARGET_CRS)
    
    # Clip Natural Earth layers
    lakes_final = gpd.clip(lakes_gdf.to_crs(TARGET_CRS), viewport_gdf)
    rivers_proj = rivers_gdf.to_crs(TARGET_CRS)
    rivers_clipped = gpd.clip(rivers_proj, viewport_gdf)
    ocean_final = gpd.clip(ocean_gdf.to_crs(TARGET_CRS), viewport_gdf)
    
    # Remove rivers that will be drawn from OSM (to avoid duplication)
    osm_names = "Péribonka|Mistassini|Ashuapmushuan|Métabetchouane|Saguenay|Ha! Ha!|Belle Rivière|Décharge"
    rivers_clean = rivers_clipped[
        ~rivers_clipped['name'].astype(str).str.contains(
            osm_names, case=False, na=False, regex=True
        )
    ]
    
    # Process OSM polygons
    water_polygons: List[Polygon] = []
    for poly in osm_polygons:
        try:
            g = gpd.GeoDataFrame(
                {'geometry': [poly]}, 
                crs=SOURCE_CRS
            ).to_crs(TARGET_CRS).iloc[0].geometry
            clipped = g.intersection(viewport_box)
            if not clipped.is_empty:
                water_polygons.append(clipped)
        except (IndexError, AttributeError):
            continue
    
    # Process OSM lines
    if osm_lines:
        all_lines = [line for lines in osm_lines.values() for line in lines]
        rivers_osm = gpd.clip(
            gpd.GeoDataFrame({'geometry': all_lines}, crs=SOURCE_CRS).to_crs(TARGET_CRS),
            viewport_gdf
        )
        rivers_final = pd.concat([rivers_clean, rivers_osm], ignore_index=True)
    else:
        rivers_final = rivers_clean
    
    return (gdf_proj, viewport_box, viewport_gdf, lakes_final, 
            rivers_final, rivers_proj, ocean_final, water_polygons)


def consolidate_river_geometries(
    osm_lines: Dict[str, List[LineString]], 
    rivers_ne: gpd.GeoDataFrame, 
    viewport_box: Polygon
) -> Dict[str, Dict[str, Union[float, Union[LineString, MultiLineString]]]]:
    """
    Consolidate river geometries from multiple sources for labeling.
    
    Parameters
    ----------
    osm_lines : dict
        OSM river centerlines by name
    rivers_ne : gpd.GeoDataFrame
        Natural Earth rivers (projected to TARGET_CRS)
    viewport_box : Polygon
        Viewport bounding box for clipping
    
    Returns
    -------
    dict
        Mapping of river names to {'total_length': float, 'geometry': LineString/MultiLineString}
    """
    river_groups: Dict[str, Dict[str, Union[float, Union[LineString, MultiLineString]]]] = {}
    
    # Process OSM lines
    for name, lines in osm_lines.items():
        proj_list: List[LineString] = []
        for line in lines:
            try:
                g = gpd.GeoDataFrame(
                    {'geometry': [line]}, 
                    crs=SOURCE_CRS
                ).to_crs(TARGET_CRS).iloc[0].geometry
                clipped = g.intersection(viewport_box)
                if not clipped.is_empty:
                    proj_list.extend(
                        clipped.geoms if isinstance(clipped, MultiLineString) else [clipped]
                    )
            except (IndexError, AttributeError):
                continue
        
        if proj_list:
            river_groups[name] = {
                'total_length': sum(l.length for l in proj_list),
                'geometry': linemerge(proj_list)
            }
    
    # Process Natural Earth rivers
    for _, row in rivers_ne.iterrows():
        name = row.get('name', '')
        if pd.isna(name):
            continue
        
        clipped = row.geometry.intersection(viewport_box)
        if clipped.is_empty:
            continue
        
        if name not in river_groups:
            river_groups[name] = {
                'total_length': 0, 
                'geometry': MultiLineString([])
            }
        
        river_groups[name]['total_length'] += clipped.length
        
        new_geoms = [clipped] if isinstance(clipped, LineString) else list(clipped.geoms)
        ex = river_groups[name]['geometry']
        if isinstance(ex, LineString):
            new_geoms.append(ex)
        elif isinstance(ex, MultiLineString):
            new_geoms.extend(ex.geoms)
        
        river_groups[name]['geometry'] = linemerge(new_geoms)
    
    return river_groups


# =============================================================================
# CIRCLE DISPERSION ALGORITHM
# =============================================================================

def disperse_points(
    x_coords: np.ndarray, 
    y_coords: np.ndarray, 
    radii: np.ndarray,
    iterations: int = DISPERSION_ITERATIONS
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Adjust point positions to prevent circle overlap using force-directed layout.
    
    This algorithm iteratively applies repulsion forces between overlapping circles
    and attraction forces toward original positions to minimize displacement while
    eliminating overlaps.
    
    Parameters
    ----------
    x_coords, y_coords : np.ndarray
        Initial coordinates of circle centers
    radii : np.ndarray
        Circle radii
    iterations : int, optional
        Number of force calculation iterations
    
    Returns
    -------
    tuple of np.ndarray
        (adjusted_x, adjusted_y) coordinates after dispersion
    """
    n = len(x_coords)
    curr_x = np.array(x_coords, dtype=float)
    curr_y = np.array(y_coords, dtype=float)
    orig_x, orig_y = curr_x.copy(), curr_y.copy()
    r_padded = np.array(radii, dtype=float) * PADDING_FACTOR
    
    for _ in range(iterations):
        fx, fy = np.zeros(n), np.zeros(n)
        
        # Calculate pairwise repulsion forces
        for i in range(n):
            for j in range(i + 1, n):
                dx, dy = curr_x[i] - curr_x[j], curr_y[i] - curr_y[j]
                dist = np.sqrt(dx * dx + dy * dy)
                min_dist = r_padded[i] + r_padded[j]
                
                if dist < min_dist:
                    overlap = min_dist - dist
                    # Avoid division by zero with random perturbation
                    nx, ny = (
                        (np.random.rand(), np.random.rand()) if dist == 0 
                        else (dx / dist, dy / dist)
                    )
                    push = overlap * REPULSION_FACTOR
                    fx[i] += nx * push
                    fy[i] += ny * push
                    fx[j] -= nx * push
                    fy[j] -= ny * push
        
        # Add attraction to original positions
        fx -= (curr_x - orig_x) * ATTRACTION_FACTOR
        fy -= (curr_y - orig_y) * ATTRACTION_FACTOR
        
        # Update positions
        curr_x += fx
        curr_y += fy
    
    return curr_x, curr_y


# =============================================================================
# TEXT LABELING UTILITIES
# =============================================================================

def get_longest_segment(geometry: Union[LineString, MultiLineString]) -> Optional[LineString]:
    """
    Extract longest continuous segment from line geometry.
    
    Parameters
    ----------
    geometry : LineString or MultiLineString
        Input line geometry
    
    Returns
    -------
    LineString or None
        Longest continuous segment, or None if geometry is empty
    """
    if isinstance(geometry, LineString):
        return geometry
    if isinstance(geometry, MultiLineString) and geometry.geoms:
        return max(geometry.geoms, key=lambda x: x.length)
    return None


def extract_line_segment(
    line: LineString, 
    start_frac: float, 
    end_frac: float
) -> Optional[LineString]:
    """
    Extract segment of line by fractional position along its length.
    
    Parameters
    ----------
    line : LineString
        Input line geometry
    start_frac, end_frac : float
        Start and end positions as fractions (0-1) of total length
    
    Returns
    -------
    LineString or None
        Extracted segment with 100 interpolated points, or None if invalid
    """
    if not line or line.length == 0:
        return None
    
    distances = np.linspace(line.length * start_frac, line.length * end_frac, 100)
    points = [line.interpolate(d) for d in distances]
    return LineString([(p.x, p.y) for p in points])


def smooth_angles(angles: np.ndarray, window: int = 3) -> np.ndarray:
    """
    Smooth angles using circular moving average.
    
    This accounts for circular nature of angles by smoothing sin/cos components
    separately to avoid discontinuities at angle wrapping points.
    
    Parameters
    ----------
    angles : np.ndarray
        Input angles in degrees
    window : int, optional
        Smoothing window size
    
    Returns
    -------
    np.ndarray
        Smoothed angles in degrees
    """
    if window < 2:
        return angles
    
    rads = np.radians(angles)
    kernel = np.ones(window) / window
    smooth_sin = np.convolve(np.sin(rads), kernel, mode='same')
    smooth_cos = np.convolve(np.cos(rads), kernel, mode='same')
    return np.degrees(np.arctan2(smooth_sin, smooth_cos))


def plot_text_along_line(
    ax: plt.Axes, 
    line: LineString, 
    text: str,
    fontsize: int = 10, 
    color: str = RIVER_LABEL_COLOR, 
    offset: float = 0,
    geom_smooth: int = 5, 
    angle_smooth: int = 5
) -> None:
    """
    Plot text following the path of a line geometry.
    
    Characters are positioned along a smoothed spline interpolation of the line,
    with orientation matching the local tangent angle. This creates curved text
    labels that follow river paths naturally.
    
    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axes to plot on
    line : LineString
        Path geometry to follow
    text : str
        Text to display
    fontsize : int, optional
        Font size in points
    color : str, optional
        Text color
    offset : float, optional
        Perpendicular offset from line in map units (positive = left of direction)
    geom_smooth : int, optional
        Window size for geometry smoothing
    angle_smooth : int, optional
        Window size for angle smoothing
    """
    if not line or line.length == 0:
        return
    
    coords = np.array(line.coords)
    # Ensure left-to-right orientation for readability
    if coords[0][0] > coords[-1][0]:
        coords = coords[::-1]
    
    # Smooth geometry to reduce noise in angle calculations
    if len(coords) > geom_smooth:
        kernel = np.ones(geom_smooth) / geom_smooth
        coords = np.column_stack([
            np.convolve(coords[:, 0], kernel, mode='valid'),
            np.convolve(coords[:, 1], kernel, mode='valid')
        ])
    
    # Calculate cumulative distances along path
    dist = np.insert(
        np.cumsum(np.sqrt(np.sum(np.diff(coords, axis=0)**2, axis=1))), 
        0, 0
    )
    if dist[-1] == 0:
        return
    
    # Create spline interpolation for smooth character placement
    try:
        k = min(3, len(coords) - 1)
        spline_x = make_interp_spline(dist, coords[:, 0], k=k)
        spline_y = make_interp_spline(dist, coords[:, 1], k=k)
    except (ValueError, TypeError):
        return
    
    # Calculate character positions along centerline
    n = len(text)
    spacing = (dist[-1] * TEXT_PATH_COVERAGE) / max(n - 1, 1)
    start = (dist[-1] - spacing * (n - 1)) / 2
    
    positions: List[Tuple[float, float]] = []
    angles: List[float] = []
    
    for i in range(n):
        d = np.clip(start + i * spacing, 0, dist[-1])
        positions.append((float(spline_x(d)), float(spline_y(d))))
        
        # Calculate tangent angle using forward/backward differences
        dp = min(dist[-1], d + ANGLE_CALC_DELTA)
        dm = max(0, d - ANGLE_CALC_DELTA)
        dx = float(spline_x(dp) - spline_x(dm))
        dy = float(spline_y(dp) - spline_y(dm))
        angles.append(math.degrees(math.atan2(dy, dx)))
    
    angles_smoothed = smooth_angles(np.array(angles), angle_smooth)
    
    # Plot individual characters
    for i, char in enumerate(text):
        x, y = positions[i]
        angle = angles_smoothed[i]
        
        if offset:
            perp = math.radians(angle + 90)
            x += offset * math.cos(perp)
            y += offset * math.sin(perp)
        
        ax.text(
            x, y, char, 
            fontsize=fontsize, 
            color=color, 
            rotation=angle,
            rotation_mode='anchor', 
            ha='center', 
            va='bottom',
            fontweight='bold', 
            fontstyle='italic', 
            zorder=7, 
            alpha=0.85
        )


def plot_horizontal_text(
    ax: plt.Axes, 
    x: float, 
    y: float, 
    text: str,
    fontsize: int = 12, 
    color: str = '#444444'
) -> None:
    """
    Plot horizontal text with white outline for visibility.
    
    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axes to plot on
    x, y : float
        Position coordinates in data units
    text : str
        Text to display
    fontsize : int, optional
        Font size in points
    color : str, optional
        Text color
    """
    ax.text(
        x, y, text, 
        fontsize=fontsize, 
        color=color, 
        ha='center', 
        va='center',
        fontweight='bold', 
        style='italic', 
        zorder=7,
        path_effects=[pe.withStroke(linewidth=3, foreground="white")]
    )


# =============================================================================
# VISUALIZATION - MAP CREATION
# =============================================================================

def create_main_map(
    gdf_proj: gpd.GeoDataFrame, 
    lakes: gpd.GeoDataFrame,
    rivers: gpd.GeoDataFrame, 
    ocean: gpd.GeoDataFrame,
    water_polygons: List[Polygon], 
    bounds: Tuple[float, float, float, float]
) -> Tuple[plt.Figure, plt.Axes]:
    """
    Create the main map figure with water features and population circles.
    
    Parameters
    ----------
    gdf_proj : gpd.GeoDataFrame
        Projected population data with 'new_x', 'new_y', 'radius_m', 'Color' columns
    lakes, rivers, ocean : gpd.GeoDataFrame
        Water feature layers (projected)
    water_polygons : list
        Additional water body polygons from OSM
    bounds : tuple
        Map extent (minx, miny, maxx, maxy) in projected coordinates
    
    Returns
    -------
    tuple
        (figure, axes) matplotlib objects
    """
    fig, ax = plt.subplots(figsize=FIGURE_SIZE)
    ax.set_facecolor(BACKGROUND_COLOR)
    ax.set_aspect('equal')
    
    # Configure frame
    for spine in ax.spines.values():
        spine.set_edgecolor('black')
        spine.set_linewidth(1.0)
    
    # Plot water features
    ocean.plot(ax=ax, color=WATER_COLOR, zorder=1)
    lakes.plot(ax=ax, color=WATER_COLOR, zorder=1)
    for poly in water_polygons:
        gpd.GeoDataFrame({'geometry': [poly]}, crs=TARGET_CRS).plot(
            ax=ax, color=WATER_COLOR, zorder=1
        )
    if not rivers.empty:
        rivers.plot(ax=ax, color=WATER_COLOR, linewidth=1.2, zorder=2)
    
    # Plot population circles
    for _, row in gdf_proj.iterrows():
        ax.add_patch(mpatches.Circle(
            (row['new_x'], row['new_y']), 
            row['radius_m'],
            facecolor=row['Color'], 
            edgecolor='black', 
            linewidth=0.5, 
            alpha=0.85, 
            zorder=4
        ))
    
    # Add city labels for larger populations
    for _, row in gdf_proj.iterrows():
        if row.get('Population size', 0) >= LABEL_SIZE_THRESHOLD:
            city = row["Location of proband's marriage"]
            if city == 'Arvida':
                continue
            
            r = row['radius_m']
            if city == 'Jonquière':
                pos = (row['new_x'], row['new_y'] - r - CITY_LABEL_OFFSET)
                align = ('center', 'top')
            else:
                pos = (row['new_x'] + r + CITY_LABEL_OFFSET, row['new_y'])
                align = ('left', 'center')
            
            ax.text(
                *pos, city, 
                ha=align[0], 
                va=align[1], 
                fontsize=8, 
                fontweight='bold',
                zorder=6, 
                path_effects=[pe.withStroke(linewidth=3, foreground="white")]
            )
    
    # Set map extent
    minx, miny, maxx, maxy = bounds
    ax.set_xlim(minx - BUFFER_METERS, maxx + BUFFER_METERS)
    ax.set_ylim(miny - BUFFER_METERS, maxy + BUFFER_METERS)
    ax.set_xticks([])
    ax.set_yticks([])
    
    # Add title
    ax.text(
        0.5, 0.02, 
        "Saguenay–Lac-Saint-Jean", 
        transform=ax.transAxes,
        ha='center', 
        va='bottom', 
        fontsize=18, 
        fontweight='bold', 
        color='black', 
        zorder=5
    )
    
    return fig, ax


def add_river_labels(
    ax: plt.Axes, 
    rivers: Dict[str, Dict[str, Union[float, Union[LineString, MultiLineString]]]], 
    inset_box: Polygon
) -> None:
    """
    Add labels to major rivers.
    
    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axes to plot on
    rivers : dict
        Consolidated river geometries from consolidate_river_geometries()
    inset_box : Polygon
        Inset map bounding box to avoid when placing labels
    """
    # Configuration: (display_text, search_pattern, start_frac, end_frac, 
    #                 fontsize, offset, avoid_inset, exclude_pattern, 
    #                 geom_smooth, angle_smooth)
    configs = [
        ("Péribonka", "Péribonka", 0.3, 0.7, 11, 1000, False, "Petite", 10, 7),
        ("Mistassini", "Mistassini", 0.05, 0.7, 11, -1000, True, None, 10, 5),
        ("Ashuapmushuan", "Ashuapmushuan", 0.00, 1.00, 10, 800, False, None, 20, 5),
        ("Métabetchouane", "Métabetchouane", 0.10, 0.90, 9, 0, False, None, 10, 3),
        ("Ha! Ha!", "Ha! Ha!", 0.1, 0.9, 9, 0, False, None, 3, 3),
        ("Belle Rivière", "Belle", 0.05, 0.95, 6, 0, False, None, 3, 3),
    ]
    
    for text, pattern, sf, ef, fs, off, avoid, excl, gs, asw in configs:
        # Find matching rivers (filter by pattern and exclusion)
        matches = [
            (n, d['total_length'], d['geometry']) 
            for n, d in rivers.items()
            if pattern.lower() in n.lower() and not (excl and excl.lower() in n.lower())
        ]
        
        if not matches:
            continue
        
        # Use longest matching river
        _, _, geom = max(matches, key=lambda x: x[1])
        path = get_longest_segment(geom)
        
        # Avoid inset box if requested
        if avoid and path:
            try:
                path = get_longest_segment(path.difference(inset_box))
            except (AttributeError, TypeError):
                pass
        
        # Only label sufficiently long segments
        if path and path.length > RIVER_LABEL_MIN_LENGTH:
            segment = extract_line_segment(path, sf, ef)
            plot_text_along_line(ax, segment, text, fs, RIVER_LABEL_COLOR, off, gs, asw)


def add_water_body_labels(
    ax: plt.Axes, 
    lakes: gpd.GeoDataFrame,
    gdf_proj: gpd.GeoDataFrame, 
    water_polygons: List[Polygon]
) -> None:
    """
    Add labels to major water bodies (Lac Saint-Jean and Saguenay Fjord).
    
    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axes to plot on
    lakes : gpd.GeoDataFrame
        Lake features
    gdf_proj : gpd.GeoDataFrame
        Projected population data (used for Saguenay positioning)
    water_polygons : list
        OSM water body polygons (used for precise Saguenay positioning)
    """
    # Label Lac Saint-Jean
    lac = lakes[lakes['name'].str.contains(
        'Saint-Jean|Lac-Saint-Jean', 
        case=False, na=False, regex=True
    )]
    if not lac.empty:
        c = lac.iloc[0].geometry.centroid
        plot_horizontal_text(
            ax, 
            c.x + LAC_SJ_LABEL_OFFSET_X, 
            c.y + LAC_SJ_LABEL_OFFSET_Y_TOP, 
            "Lac", 
            14, 
            WATER_LABEL_COLOR
        )
        plot_horizontal_text(
            ax, 
            c.x + LAC_SJ_LABEL_OFFSET_X, 
            c.y + LAC_SJ_LABEL_OFFSET_Y_BOTTOM, 
            "Saint-Jean", 
            14, 
            WATER_LABEL_COLOR
        )
    
    # Label Saguenay Fjord
    subset = gdf_proj[gdf_proj['Color'] == SUBDIVISION_COLORS["East of Ha! Ha!"]]
    if not subset.empty:
        # Position label based on "East of Ha! Ha!" subdivision
        avg_x = (subset.geometry.x.min() + subset.geometry.x.max()) / 2
        max_y = subset.geometry.y.max()
        final_x, final_y = avg_x, max_y + SAGUENAY_DEFAULT_OFFSET_Y
        
        # Try to find exact water polygon position using ray casting
        if water_polygons:
            ray = LineString([
                (avg_x, max_y), 
                (avg_x, max_y + SAGUENAY_RAY_LENGTH)
            ])
            intersection = unary_union(water_polygons).intersection(ray)
            
            if not intersection.is_empty:
                if isinstance(intersection, MultiLineString):
                    seg = min(intersection.geoms, key=lambda g: g.centroid.y)
                elif isinstance(intersection, LineString):
                    seg = intersection
                else:
                    seg = None
                
                if seg:
                    final_x, final_y = seg.centroid.x, seg.centroid.y
        
        plot_horizontal_text(
            ax, 
            final_x, 
            final_y + SAGUENAY_LABEL_OFFSET_Y_TOP, 
            "Saguenay", 
            14, 
            WATER_LABEL_COLOR
        )
        plot_horizontal_text(
            ax, 
            final_x, 
            final_y + SAGUENAY_LABEL_OFFSET_Y_BOTTOM, 
            "Fjord", 
            14, 
            WATER_LABEL_COLOR
        )


def create_inset_map(
    fig: plt.Figure, 
    ax: plt.Axes, 
    admin: gpd.GeoDataFrame,
    lakes: gpd.GeoDataFrame, 
    ocean: gpd.GeoDataFrame,
    rivers: gpd.GeoDataFrame, 
    viewport: gpd.GeoDataFrame
) -> None:
    """
    Create Quebec context inset map showing main map extent.
    
    Parameters
    ----------
    fig : matplotlib.figure.Figure
        Main figure
    ax : matplotlib.axes.Axes
        Main axes to attach inset to
    admin : gpd.GeoDataFrame
        Administrative boundaries (projected)
    lakes, ocean, rivers : gpd.GeoDataFrame
        Water features (projected)
    viewport : gpd.GeoDataFrame
        Main map viewport to highlight
    """
    ax_ins = inset_axes(
        ax, 
        width=f"{INSET_SIZE_FRACTION * 100}%", 
        height=f"{INSET_SIZE_FRACTION * 100}%", 
        loc='upper left', 
        borderpad=1
    )
    ax_ins.set_facecolor(WATER_COLOR)
    
    for spine in ax_ins.spines.values():
        spine.set_edgecolor('black')
        spine.set_linewidth(1.0)
    
    # Extract Quebec boundaries
    admin_proj = admin.to_crs(TARGET_CRS)
    canada = admin_proj[
        admin_proj.get('admin', admin_proj.get('adm0_name')) == 'Canada'
    ]
    quebec = canada[canada['name'].astype(str).str.contains(
        r'Qu.bec', regex=True, case=False
    )]
    if quebec.empty:
        quebec = canada[canada['name'] == 'Quebec']
    
    # Calculate inset extent with padding
    minx, miny, maxx, maxy = quebec.total_bounds
    pad_x = (maxx - minx) * INSET_PADDING_FRACTION
    pad_y = (maxy - miny) * INSET_PADDING_FRACTION
    view_poly = box(minx - pad_x, miny - pad_y, maxx + pad_x, maxy + pad_y)
    view_gdf = gpd.GeoDataFrame({'geometry': [view_poly]}, crs=TARGET_CRS)
    
    # Plot administrative boundaries
    gpd.clip(
        admin_proj[~admin_proj.index.isin(quebec.index)], 
        view_gdf
    ).plot(
        ax=ax_ins, 
        color='#d3d3d3', 
        edgecolor='white', 
        linewidth=0.4, 
        zorder=1
    )
    quebec.plot(
        ax=ax_ins, 
        color='white', 
        edgecolor='#666666', 
        linewidth=0.4, 
        zorder=1
    )
    
    # Add water features
    gpd.clip(lakes, view_gdf).plot(ax=ax_ins, color=WATER_COLOR, zorder=2)
    gpd.clip(ocean, view_gdf).plot(ax=ax_ins, color=WATER_COLOR, zorder=2)
    gpd.clip(rivers, view_gdf).plot(
        ax=ax_ins, 
        color=WATER_COLOR, 
        linewidth=0.6, 
        zorder=2
    )
    
    # Highlight main map extent
    viewport.plot(
        ax=ax_ins, 
        facecolor='none', 
        edgecolor='red', 
        linewidth=1.0, 
        zorder=3
    )
    
    ax_ins.set_xlim(minx - pad_x, maxx + pad_x)
    ax_ins.set_ylim(miny - pad_y, maxy + pad_y)
    ax_ins.set_xticks([])
    ax_ins.set_yticks([])
    ax_ins.text(
        0.5, 0.05, 
        "Québec", 
        transform=ax_ins.transAxes,
        fontsize=16, 
        fontweight='bold', 
        ha='center', 
        va='bottom', 
        zorder=4
    )


def add_legends(fig: plt.Figure, ax: plt.Axes, gdf_proj: gpd.GeoDataFrame) -> None:
    """
    Add color-coded subdivision legend and population size legend.
    
    Parameters
    ----------
    fig : matplotlib.figure.Figure
        Figure object
    ax : matplotlib.axes.Axes
        Main axes
    gdf_proj : gpd.GeoDataFrame
        Projected population data (used for coordinate transformation)
    """
    # Color legend for subdivisions
    patches = [
        mpatches.Patch(color=c, label=l) 
        for l, c in SUBDIVISION_COLORS.items()
    ]
    leg1 = ax.legend(
        handles=patches, 
        loc='upper right', 
        title="Watercourse Subdivisions",
        fontsize=10, 
        framealpha=1.0, 
        bbox_to_anchor=(0.98, 0.98),
        edgecolor='black', 
        title_fontproperties={'weight': 'bold', 'size': 12}
    )
    ax.add_artist(leg1)
    
    # Calculate position for size legend below color legend
    plt.tight_layout()
    fig.canvas.draw()
    
    trans_inv = fig.transFigure.inverted()
    bbox = leg1.get_window_extent().transformed(trans_inv)
    
    # Convert map units (meters) to figure units for accurate circle sizing
    # We measure the figure-space distance corresponding to 10km in data space
    p0 = ax.transData.transform((0, 0))
    p1 = ax.transData.transform((10000, 0))
    fig_per_m = (
        (trans_inv.transform((p1[0] - p0[0], 0))[0] - trans_inv.transform((0, 0))[0]) 
        / 10000
    )
    
    # Population size legend
    radii = [
        np.sqrt(p) * POPULATION_SCALE_FACTOR * fig_per_m 
        for p in LEGEND_POPULATION_SIZES
    ]
    
    leg2_h = (LEGEND_PAD + LEGEND_TITLE_HEIGHT + LEGEND_SPACING + 
              max(radii) * 2 + LEGEND_SPACING + LEGEND_LABEL_HEIGHT + LEGEND_PAD)
    leg2_x, leg2_y = bbox.x0, bbox.y0 - LEGEND_GAP - leg2_h
    
    # Background
    ax.add_patch(mpatches.FancyBboxPatch(
        (leg2_x, leg2_y), 
        bbox.width, 
        leg2_h,
        boxstyle="round,pad=0,rounding_size=0.01",
        facecolor='white', 
        edgecolor='black', 
        linewidth=1, 
        alpha=1.0,
        zorder=20, 
        transform=fig.transFigure
    ))
    
    # Title
    ax.text(
        leg2_x + bbox.width / 2, 
        leg2_y + leg2_h - LEGEND_PAD, 
        "Population Size",
        ha='center', 
        va='top', 
        fontsize=12, 
        fontweight='bold',
        transform=fig.transFigure, 
        zorder=21
    )
    
    # Circles and labels
    col_w = (bbox.width - 2 * LEGEND_PAD) / 3
    circle_y = leg2_y + LEGEND_PAD + LEGEND_LABEL_HEIGHT + LEGEND_SPACING
    
    for i, (label, r) in enumerate(zip(LEGEND_POPULATION_LABELS, radii)):
        cx = leg2_x + LEGEND_PAD + i * col_w + col_w / 2
        ax.add_patch(mpatches.Circle(
            (cx, circle_y + r), 
            r, 
            facecolor='gray',
            edgecolor='black', 
            linewidth=0.5, 
            alpha=0.7,
            transform=fig.transFigure, 
            zorder=21
        ))
        ax.text(
            cx, 
            leg2_y + LEGEND_PAD + LEGEND_LABEL_HEIGHT / 2, 
            label, 
            ha='center', 
            va='center',
            fontsize=10, 
            transform=fig.transFigure, 
            zorder=21
        )


# =============================================================================
# MAIN EXECUTION
# =============================================================================

def main() -> None:
    """
    Main execution function.
    
    Orchestrates the complete map generation pipeline:
    1. Load population data
    2. Download geographic base layers
    3. Query OpenStreetMap for detailed river data
    4. Project and clip all layers to study area
    5. Consolidate river geometries from multiple sources
    6. Calculate circle sizes and adjust positions to prevent overlap
    7. Create main map with water features and population circles
    8. Add river and water body labels
    9. Create Quebec context inset map
    10. Add legends
    11. Save output as SVG and PNG
    """
    # Load data
    gdf = load_population_data("CharlevoixGradient.csv")
    admin, lakes, rivers_ne, ocean = download_natural_earth_data()
    
    # Calculate bounding box for OSM query
    minx, miny, maxx, maxy = gdf.total_bounds
    bbox = (
        miny - BUFFER_DEGREES, minx - BUFFER_DEGREES,
        maxy + BUFFER_DEGREES, maxx + BUFFER_DEGREES
    )
    
    # Query OpenStreetMap
    osm_lines, osm_polygons = query_openstreetmap_rivers(bbox)
    
    # Project and clip layers
    (gdf_proj, viewport_box, viewport_gdf, lakes_final, rivers_final,
     rivers_proj, ocean_final, water_polygons) = project_and_clip_layers(
        gdf, lakes, rivers_ne, ocean, osm_lines, osm_polygons
    )
    
    # Consolidate river geometries for labeling
    river_groups = consolidate_river_geometries(osm_lines, rivers_proj, viewport_box)
    
    # Calculate circle sizes and adjust positions
    print("Adjusting circle positions...")
    gdf_proj['radius_m'] = np.sqrt(gdf_proj['Population size']) * POPULATION_SCALE_FACTOR
    gdf_proj['new_x'], gdf_proj['new_y'] = disperse_points(
        gdf_proj.geometry.x.values, 
        gdf_proj.geometry.y.values, 
        gdf_proj['radius_m'].values
    )
    
    # Create main map
    fig, ax = create_main_map(
        gdf_proj, lakes_final, rivers_final, ocean_final,
        water_polygons, gdf_proj.total_bounds
    )
    
    # Calculate inset box for river label placement
    minx, miny, maxx, maxy = gdf_proj.total_bounds
    inset_box = box(
        minx - BUFFER_METERS,
        maxy + BUFFER_METERS - (maxy - miny + 2 * BUFFER_METERS) * INSET_SIZE_FRACTION,
        minx - BUFFER_METERS + (maxx - minx + 2 * BUFFER_METERS) * INSET_SIZE_FRACTION,
        maxy + BUFFER_METERS
    )
    
    # Add labels
    add_river_labels(ax, river_groups, inset_box)
    add_water_body_labels(ax, lakes_final, gdf_proj, water_polygons)
    
    # Create inset map
    create_inset_map(
        fig, ax, admin, 
        lakes.to_crs(TARGET_CRS),
        ocean.to_crs(TARGET_CRS), 
        rivers_ne.to_crs(TARGET_CRS), 
        viewport_gdf
    )
    
    # Add legends
    add_legends(fig, ax, gdf_proj)
    
    # Save outputs
    plt.savefig("map_slsj.svg", dpi=300, format='svg', bbox_inches='tight')
    plt.savefig("map_slsj.png", dpi=300, format='png', bbox_inches='tight')
    print("Saved map_slsj.svg and .png")
    plt.show()


if __name__ == "__main__":
    main()