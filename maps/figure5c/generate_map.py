"""
Saguenay–Lac-Saint-Jean Genetic Contribution Map Generator
Complete version with proper OSM filtering, dispersion, and labeling

Author: Derived from the Figure 1 map code using Claude 4.5 Sonnet
"""

import warnings
import zipfile
import tempfile
import math
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Union

import numpy as np
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.patheffects as pe
from shapely.geometry import box, LineString, MultiLineString, Polygon
from shapely.ops import linemerge, unary_union
from matplotlib.patches import Wedge
from scipy.interpolate import make_interp_spline

warnings.filterwarnings('ignore')

# =============================================================================
# CONFIGURATION
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

CHARLEVOIX_COLORS: Dict[str, str] = {
    "La Malbaie": 'tab:blue',
    "Baie-Saint-Paul": 'tab:orange',
    "Les Éboulements": 'tab:green',
    "Elsewhere in Charlevoix": 'tab:red',
    "Not in Charlevoix": 'tab:gray'
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

# Visualization parameters
TARGET_CRS = "EPSG:32198"
SOURCE_CRS = "EPSG:4326"
WATER_COLOR = '#a3ccff'
WATER_LABEL_COLOR = '#3a6ea5'
RIVER_LABEL_COLOR = '#555555'
BACKGROUND_COLOR = 'white'
BUFFER_METERS = 25000
POPULATION_SCALE_FACTOR = 125.0  # 5x larger than 25
FIGURE_SIZE = (10, 10)

# Dispersion parameters
DISPERSION_ITERATIONS = 50
REPULSION_FACTOR = 0.5
ATTRACTION_FACTOR = 0.05
PADDING_FACTOR = 1.05

# Manual position adjustments (in meters) - adjust these to fine-tune circle positions
# Format: subdivision_name: (x_offset, y_offset)
# Positive x = eastward, negative x = westward
# Positive y = northward, negative y = southward
MANUAL_POSITION_ADJUSTMENTS: Dict[str, Tuple[float, float]] = {
    "East of Ha! Ha!": (0, 0),
    "South of Saguenay River": (-15000, 0),
    "North of Saguenay River": (0, 0),
    "East of Lac Saint-Jean": (30000, -10000),
    "La-Belle-Rivière–Métabetchouane": (5000, 7000),
    "Mistassini-Péribonka": (-12000, 0),
    "Métabetchouane-Ashuapmushuan": (19000, -16000),
    "Ashuapmushuan-Mistassini": (22000, 4000),
}

# Label offsets
LAC_SJ_LABEL_OFFSET_X = 3000
LAC_SJ_LABEL_OFFSET_Y_TOP = 5000
LAC_SJ_LABEL_OFFSET_Y_BOTTOM = -3000
SAGUENAY_LABEL_OFFSET_Y_TOP = 5000
SAGUENAY_LABEL_OFFSET_Y_BOTTOM = -3000
SAGUENAY_DEFAULT_OFFSET_Y = 6000
SAGUENAY_RAY_LENGTH = 20000
RIVER_LABEL_MIN_LENGTH = 2000

# Legend parameters
LEGEND_GAP = 0.015
LEGEND_PAD = 0.015
LEGEND_TITLE_HEIGHT = 0.025
LEGEND_LABEL_HEIGHT = 0.020
LEGEND_SPACING = 0.010

# =============================================================================
# DATA LOADING
# =============================================================================

def load_geodata_from_zip(zip_path: str) -> Tuple:
    """Load geographic data from zip file."""
    print("Loading geographic data from zip...")
    
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        
        with zipfile.ZipFile(zip_path, 'r') as zip_ref:
            zip_ref.extractall(tmpdir)
        
        admin = gpd.read_file(tmpdir / 'admin')
        lakes = gpd.read_file(tmpdir / 'lakes')
        rivers = gpd.read_file(tmpdir / 'rivers')
        ocean = gpd.read_file(tmpdir / 'ocean')
        
        print(f"  ✓ Loaded admin: {len(admin)} features")
        print(f"  ✓ Loaded lakes: {len(lakes)} features")
        print(f"  ✓ Loaded rivers: {len(rivers)} features")
        print(f"  ✓ Loaded ocean: {len(ocean)} features")
        
        osm_rivers = None
        osm_water = None
        
        if (tmpdir / 'osm_rivers').exists():
            osm_rivers = gpd.read_file(tmpdir / 'osm_rivers')
            print(f"  ✓ Loaded OSM rivers: {len(osm_rivers)} features")
        
        if (tmpdir / 'osm_water').exists():
            osm_water = gpd.read_file(tmpdir / 'osm_water')
            print(f"  ✓ Loaded OSM water: {len(osm_water)} features")
        
        return admin, lakes, rivers, ocean, osm_rivers, osm_water


def load_population_data(filepath: str) -> gpd.GeoDataFrame:
    """Load population data."""
    df = pd.read_csv(filepath)
    df['Subdivision'] = df["Location of proband's marriage"].map(CITY_SUBDIVISIONS)
    df['Population size'] = df.get('Population size', 1000).fillna(500)
    
    return gpd.GeoDataFrame(
        df,
        geometry=gpd.points_from_xy(df.Longitude, df.Latitude),
        crs=SOURCE_CRS
    )


def load_genetic_contributions(filepath: str) -> pd.DataFrame:
    """Load genetic contribution data."""
    df = pd.read_csv(filepath)
    df = df.set_index('Watercourse')
    return df


def aggregate_by_subdivision(gdf: gpd.GeoDataFrame, gdf_proj: gpd.GeoDataFrame,
                             genetic_df: pd.DataFrame) -> pd.DataFrame:
    """Aggregate municipalities by watercourse subdivision."""
    results = []
    
    for subdivision in SUBDIVISION_COLORS.keys():
        sub_data = gdf_proj[gdf_proj['Subdivision'] == subdivision]
        
        if len(sub_data) == 0:
            continue
        
        total_pop = sub_data['Population size'].sum()
        weights = sub_data['Population size'].values
        centroid_x = np.average(sub_data.geometry.x.values, weights=weights)
        centroid_y = np.average(sub_data.geometry.y.values, weights=weights)
        
        if subdivision not in genetic_df.index:
            continue
            
        genetic_row = genetic_df.loc[subdivision]
        
        results.append({
            'subdivision': subdivision,
            'total_pop': total_pop,
            'centroid_x': centroid_x,
            'centroid_y': centroid_y,
            'La Malbaie': genetic_row['La Malbaie'],
            'Baie-Saint-Paul': genetic_row['Baie-Saint-Paul'],
            'Les Éboulements': genetic_row['Les Éboulements'],
            'Elsewhere in Charlevoix': genetic_row['Elsewhere in Charlevoix'],
            'Not in Charlevoix': genetic_row['Not in Charlevoix']
        })
    
    return pd.DataFrame(results)


# =============================================================================
# PROJECTION AND CLIPPING
# =============================================================================

def project_and_clip_layers(gdf: gpd.GeoDataFrame, lakes: gpd.GeoDataFrame,
                            rivers: gpd.GeoDataFrame, ocean: gpd.GeoDataFrame,
                            osm_rivers: Optional[gpd.GeoDataFrame],
                            osm_water: Optional[gpd.GeoDataFrame]) -> Tuple:
    """Project all layers to target CRS and clip to study area."""
    print("Projecting and clipping layers...")
    
    gdf_proj = gdf.to_crs(TARGET_CRS)
    
    minx, miny, maxx, maxy = gdf_proj.total_bounds
    viewport_box = box(
        minx - BUFFER_METERS, miny - BUFFER_METERS,
        maxx + BUFFER_METERS, maxy + BUFFER_METERS
    )
    
    lakes_proj = lakes.to_crs(TARGET_CRS)
    lakes_final = gpd.clip(lakes_proj, viewport_box)
    
    rivers_proj = rivers.to_crs(TARGET_CRS)
    rivers_final = gpd.clip(rivers_proj, viewport_box)
    
    ocean_proj = ocean.to_crs(TARGET_CRS)
    ocean_final = gpd.clip(ocean_proj, viewport_box)
    
    osm_rivers_final = None
    osm_water_final = []
    
    if osm_rivers is not None:
        osm_rivers_proj = osm_rivers.to_crs(TARGET_CRS)
        osm_rivers_final = gpd.clip(osm_rivers_proj, viewport_box)
        print(f"  ✓ Clipped OSM rivers: {len(osm_rivers_final)} features")
    
    if osm_water is not None:
        osm_water_proj = osm_water.to_crs(TARGET_CRS)
        osm_water_clipped = gpd.clip(osm_water_proj, viewport_box)
        osm_water_final = list(osm_water_clipped.geometry)
        print(f"  ✓ Clipped OSM water: {len(osm_water_final)} polygons")
    
    return (gdf_proj, viewport_box, lakes_final, rivers_final, rivers_proj,
            ocean_final, osm_rivers_final, osm_water_final)


# =============================================================================
# RIVER CONSOLIDATION
# =============================================================================

def consolidate_river_geometries(osm_rivers: Optional[gpd.GeoDataFrame],
                                 rivers_ne: gpd.GeoDataFrame,
                                 viewport_box: Polygon) -> Dict:
    """Consolidate river geometries for labeling."""
    river_groups = {}
    
    # Process OSM rivers
    if osm_rivers is not None:
        for _, row in osm_rivers.iterrows():
            name = row.get('name', 'Unknown')
            geom = row.geometry
            
            if name not in river_groups:
                river_groups[name] = {
                    'total_length': 0,
                    'geometry': MultiLineString([])
                }
            
            river_groups[name]['total_length'] += geom.length
            
            # Merge geometries
            existing = river_groups[name]['geometry']
            if isinstance(existing, MultiLineString) and existing.is_empty:
                river_groups[name]['geometry'] = geom
            else:
                if isinstance(geom, LineString):
                    new_geoms = [geom]
                else:
                    new_geoms = list(geom.geoms)
                
                if isinstance(existing, LineString):
                    new_geoms.append(existing)
                elif isinstance(existing, MultiLineString):
                    new_geoms.extend(existing.geoms)
                
                river_groups[name]['geometry'] = linemerge(new_geoms)
    
    return river_groups


# =============================================================================
# DISPERSION ALGORITHM
# =============================================================================

def disperse_points(x_coords: np.ndarray, y_coords: np.ndarray, 
                   radii: np.ndarray, agg_data: pd.DataFrame,
                   iterations: int = 100) -> Tuple:
    """Adjust point positions to prevent circle overlap and avoid text labels."""
    n = len(x_coords)
    curr_x = np.array(x_coords, dtype=float)
    curr_y = np.array(y_coords, dtype=float)
    orig_x, orig_y = curr_x.copy(), curr_y.copy()
    r_padded = np.array(radii, dtype=float) * PADDING_FACTOR
    
    # Calculate actual text label positions based on the data
    text_labels = []
    
    # Lac Saint-Jean labels
    lac_x = curr_x.mean()
    lac_y = curr_y.mean()
    text_labels.append((lac_x + LAC_SJ_LABEL_OFFSET_X, lac_y + LAC_SJ_LABEL_OFFSET_Y_TOP, 18000))
    text_labels.append((lac_x + LAC_SJ_LABEL_OFFSET_X, lac_y + LAC_SJ_LABEL_OFFSET_Y_BOTTOM, 18000))
    
    # Saguenay Fjord labels (based on East of Ha! Ha! position)
    east_haha = agg_data[agg_data['subdivision'] == "East of Ha! Ha!"]
    if not east_haha.empty:
        avg_x = (east_haha['centroid_x'].min() + east_haha['centroid_x'].max()) / 2
        max_y = east_haha['centroid_y'].max()
        sag_x = avg_x
        sag_y = max_y + SAGUENAY_DEFAULT_OFFSET_Y
        text_labels.append((sag_x, sag_y + SAGUENAY_LABEL_OFFSET_Y_TOP, 18000))
        text_labels.append((sag_x, sag_y + SAGUENAY_LABEL_OFFSET_Y_BOTTOM, 18000))
    
    # River label areas - more targeted positions
    # Use actual subdivision centroids to estimate river positions
    for _, row in agg_data.iterrows():
        sub = row['subdivision']
        cx, cy = row['centroid_x'], row['centroid_y']
        
        # Add repulsion zones near specific subdivisions where river labels appear
        if sub == "Mistassini-Péribonka":
            text_labels.append((cx + 15000, cy + 20000, 22000))  # Péribonka label area
            text_labels.append((cx - 10000, cy + 10000, 20000))  # Mistassini label area
        elif sub == "Ashuapmushuan-Mistassini":
            text_labels.append((cx + 10000, cy, 25000))  # Ashuapmushuan label area
        elif sub == "Métabetchouane-Ashuapmushuan":
            text_labels.append((cx + 10000, cy - 8000, 18000))  # Métabetchouane label area
        elif sub == "East of Ha! Ha!":
            text_labels.append((cx, cy - 10000, 15000))  # Ha! Ha! label area
        elif sub == "La-Belle-Rivière–Métabetchouane":
            text_labels.append((cx + 15000, cy, 12000))  # Belle Rivière label area
    
    for iteration in range(iterations):
        fx, fy = np.zeros(n), np.zeros(n)
        
        # Calculate pairwise repulsion forces between circles
        for i in range(n):
            for j in range(i + 1, n):
                dx, dy = curr_x[i] - curr_x[j], curr_y[i] - curr_y[j]
                dist = np.sqrt(dx * dx + dy * dy)
                min_dist = r_padded[i] + r_padded[j]
                
                if dist < min_dist:
                    overlap = min_dist - dist
                    nx, ny = (
                        (np.random.rand(), np.random.rand()) if dist == 0 
                        else (dx / dist, dy / dist)
                    )
                    push = overlap * REPULSION_FACTOR
                    fx[i] += nx * push
                    fy[i] += ny * push
                    fx[j] -= nx * push
                    fy[j] -= ny * push
        
        # Add repulsion from text label positions
        label_repulsion_factor = 0.4  # Stronger repulsion from labels
        for i in range(n):
            for label_x, label_y, label_radius in text_labels:
                dx = curr_x[i] - label_x
                dy = curr_y[i] - label_y
                dist = np.sqrt(dx * dx + dy * dy)
                min_dist = r_padded[i] + label_radius
                
                if dist < min_dist:
                    overlap = min_dist - dist
                    nx, ny = (
                        (np.random.rand(), np.random.rand()) if dist == 0 
                        else (dx / dist, dy / dist)
                    )
                    push = overlap * label_repulsion_factor
                    fx[i] += nx * push
                    fy[i] += ny * push
        
        # Add attraction to original positions (weaker to allow more movement)
        fx -= (curr_x - orig_x) * (ATTRACTION_FACTOR * 0.5)
        fy -= (curr_y - orig_y) * (ATTRACTION_FACTOR * 0.5)
        
        # Update positions
        curr_x += fx
        curr_y += fy
    
    return curr_x, curr_y


# =============================================================================
# TEXT LABELING UTILITIES
# =============================================================================

def get_longest_segment(geometry: Union[LineString, MultiLineString]) -> Optional[LineString]:
    """Extract longest continuous segment from line geometry."""
    if isinstance(geometry, LineString):
        return geometry
    if isinstance(geometry, MultiLineString) and geometry.geoms:
        return max(geometry.geoms, key=lambda x: x.length)
    return None


def extract_line_segment(line: LineString, start_frac: float, end_frac: float) -> Optional[LineString]:
    """Extract segment of line by fractional position along its length."""
    if not line or line.length == 0:
        return None
    
    distances = np.linspace(line.length * start_frac, line.length * end_frac, 100)
    points = [line.interpolate(d) for d in distances]
    return LineString([(p.x, p.y) for p in points])


def smooth_angles(angles: np.ndarray, window: int = 3) -> np.ndarray:
    """Smooth angles using circular moving average."""
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
    """Plot text following the path of a line geometry (zigzag style)."""
    if not line or line.length == 0:
        return
    
    coords = np.array(line.coords)
    # Ensure left-to-right orientation
    if coords[0][0] > coords[-1][0]:
        coords = coords[::-1]
    
    # Smooth geometry
    if len(coords) > geom_smooth:
        kernel = np.ones(geom_smooth) / geom_smooth
        coords = np.column_stack([
            np.convolve(coords[:, 0], kernel, mode='valid'),
            np.convolve(coords[:, 1], kernel, mode='valid')
        ])
    
    # Calculate cumulative distances
    dist = np.insert(
        np.cumsum(np.sqrt(np.sum(np.diff(coords, axis=0)**2, axis=1))), 
        0, 0
    )
    if dist[-1] == 0:
        return
    
    # Create spline interpolation
    try:
        from scipy.interpolate import make_interp_spline
        k = min(3, len(coords) - 1)
        spline_x = make_interp_spline(dist, coords[:, 0], k=k)
        spline_y = make_interp_spline(dist, coords[:, 1], k=k)
    except (ValueError, TypeError, ImportError):
        return
    
    # Calculate character positions
    n = len(text)
    TEXT_PATH_COVERAGE = 0.8
    spacing = (dist[-1] * TEXT_PATH_COVERAGE) / max(n - 1, 1)
    start = (dist[-1] - spacing * (n - 1)) / 2
    
    positions: List[Tuple[float, float]] = []
    angles: List[float] = []
    
    ANGLE_CALC_DELTA = 100
    for i in range(n):
        d = np.clip(start + i * spacing, 0, dist[-1])
        positions.append((float(spline_x(d)), float(spline_y(d))))
        
        # Calculate tangent angle
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
    """Plot horizontal text with white outline for visibility."""
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

def draw_pie_chart(ax: plt.Axes, center_x: float, center_y: float, radius: float,
                   sizes: List[float], colors: List[str], edge_color: str,
                   start_angle: float = 90) -> None:
    """Draw a pie chart as individual wedges using actual percentages."""
    total = sum(sizes)
    if total == 0:
        return
    
    # Use actual percentages, not normalized
    current_angle = start_angle
    
    for size, color in zip(sizes, colors):
        if size <= 0:
            continue
        
        angle_extent = (size / 100) * 360
        
        wedge = Wedge(
            center=(center_x, center_y),
            r=radius,
            theta1=current_angle,
            theta2=current_angle + angle_extent,
            facecolor=color,
            edgecolor=edge_color,
            linewidth=2,
            alpha=0.9,
            zorder=10
        )
        ax.add_patch(wedge)
        
        current_angle += angle_extent


# =============================================================================
# MAP CREATION
# =============================================================================

def create_main_map(agg_data: pd.DataFrame, lakes: gpd.GeoDataFrame,
                   rivers: gpd.GeoDataFrame, ocean: gpd.GeoDataFrame,
                   osm_rivers: Optional[gpd.GeoDataFrame],
                   osm_water: List[Polygon],
                   rivers_proj: gpd.GeoDataFrame,
                   viewport_box: Polygon,
                   bounds: Tuple[float, float, float, float]) -> Tuple[plt.Figure, plt.Axes]:
    """Create the main map figure."""
    fig, ax = plt.subplots(figsize=FIGURE_SIZE)
    ax.set_facecolor(BACKGROUND_COLOR)
    ax.set_aspect('equal')
    
    for spine in ax.spines.values():
        spine.set_edgecolor('black')
        spine.set_linewidth(1.0)
    
    # Plot water features
    ocean.plot(ax=ax, color=WATER_COLOR, zorder=1)
    lakes.plot(ax=ax, color=WATER_COLOR, zorder=1)
    if not rivers.empty:
        rivers.plot(ax=ax, color=WATER_COLOR, linewidth=1.2, zorder=2)
    
    # Plot OSM water polygons efficiently
    if osm_water:
        print(f"  Plotting {len(osm_water)} OSM water polygons...")
        osm_water_gdf = gpd.GeoDataFrame({'geometry': osm_water}, crs=TARGET_CRS)
        osm_water_gdf.plot(ax=ax, color=WATER_COLOR, zorder=1)
    
    # Plot OSM rivers
    if osm_rivers is not None and not osm_rivers.empty:
        osm_rivers.plot(ax=ax, color=WATER_COLOR, linewidth=1.5, zorder=2)
    
    # Draw pie charts
    charlevoix_cols = ['La Malbaie', 'Baie-Saint-Paul', 'Les Éboulements', 
                       'Elsewhere in Charlevoix', 'Not in Charlevoix']
    
    for _, row in agg_data.iterrows():
        sizes = [row[col] for col in charlevoix_cols]
        colors = [CHARLEVOIX_COLORS[col] for col in charlevoix_cols]
        edge_color = SUBDIVISION_COLORS[row['subdivision']]
        
        radius = np.sqrt(row['total_pop']) * POPULATION_SCALE_FACTOR
        
        draw_pie_chart(
            ax,
            row['new_x'],
            row['new_y'],
            radius,
            sizes,
            colors,
            edge_color
        )
    
    # Consolidate rivers and add labels
    print("Consolidating river geometries...")
    river_groups = consolidate_river_geometries(osm_rivers, rivers_proj, viewport_box)
    
    # Create empty inset box (no inset in this version)
    inset_box = box(0, 0, 0, 0)
    add_river_labels(ax, river_groups, inset_box)
    
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


def add_river_labels(ax: plt.Axes, river_groups: Dict, inset_box: Polygon) -> None:
    """Add labels to major rivers with zigzag text style."""
    # Configuration: (display_text, search_pattern, start_frac, end_frac, 
    #                 fontsize, offset, avoid_inset, exclude_pattern, 
    #                 geom_smooth, angle_smooth)
    configs = [
        ("Péribonka", "Péribonka", 0.05, 0.95, 11, 1000, False, "Petite", 10, 7),
        ("Mistassini", "Mistassini", 0.05, 0.7, 11, -1000, True, None, 10, 5),
        ("Ashuapmushuan", "Ashuapmushuan", 0.00, 1.00, 10, 800, False, None, 20, 5),
        ("Métabetchouane", "Métabetchouane", 0.10, 0.90, 9, 0, False, None, 10, 3),
        ("Ha! Ha!", "Ha! Ha!", 0.1, 0.9, 9, 0, False, None, 3, 3),
        ("Belle Rivière", "Belle", 0.05, 0.95, 6, 0, False, None, 3, 3),
    ]
    
    for text, pattern, sf, ef, fs, off, avoid, excl, gs, asw in configs:
        # Find matching rivers
        matches = [
            (n, d['total_length'], d['geometry']) 
            for n, d in river_groups.items()
            if pattern.lower() in n.lower() and not (excl and excl.lower() in n.lower())
        ]
        
        if not matches:
            continue
        
        # Use longest matching river
        _, _, geom = max(matches, key=lambda x: x[1])
        path = get_longest_segment(geom)
        
        # Avoid inset box if requested
        if avoid and path and inset_box:
            try:
                path = get_longest_segment(path.difference(inset_box))
            except (AttributeError, TypeError):
                pass
        
        # Only label sufficiently long segments
        if path and path.length > RIVER_LABEL_MIN_LENGTH:
            segment = extract_line_segment(path, sf, ef)
            if segment:
                plot_text_along_line(ax, segment, text, fs, RIVER_LABEL_COLOR, off, gs, asw)


def add_water_labels(ax: plt.Axes, lakes: gpd.GeoDataFrame,
                    agg_data: pd.DataFrame, osm_water: List[Polygon]) -> None:
    """Add labels for Lac Saint-Jean and Saguenay Fjord."""
    # Lac Saint-Jean
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
    
    # Saguenay Fjord
    subset = agg_data[agg_data['subdivision'] == "East of Ha! Ha!"]
    if not subset.empty:
        avg_x = (subset['centroid_x'].min() + subset['centroid_x'].max()) / 2
        max_y = subset['centroid_y'].max()
        final_x, final_y = avg_x, max_y + SAGUENAY_DEFAULT_OFFSET_Y
        
        # Try to find exact position using ray casting
        if osm_water:
            ray = LineString([
                (avg_x, max_y),
                (avg_x, max_y + SAGUENAY_RAY_LENGTH)
            ])
            intersection = unary_union(osm_water).intersection(ray)
            
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


def add_legends(fig: plt.Figure, ax: plt.Axes, agg_data: pd.DataFrame) -> None:
    """Add three legends."""
    patches1 = [
        mpatches.Patch(color=c, label=l)
        for l, c in SUBDIVISION_COLORS.items()
    ]
    leg1 = ax.legend(
        handles=patches1,
        loc='upper right',
        title="Watercourse Subdivisions",
        fontsize=10,
        framealpha=1.0,
        bbox_to_anchor=(0.98, 0.98),
        edgecolor='black',
        title_fontproperties={'weight': 'bold', 'size': 12}
    )
    ax.add_artist(leg1)
    
    plt.tight_layout()
    fig.canvas.draw()
    
    trans_inv = fig.transFigure.inverted()
    bbox1 = leg1.get_window_extent().transformed(trans_inv)
    
    patches2 = [
        mpatches.Patch(color=c, label=l)
        for l, c in CHARLEVOIX_COLORS.items()
    ]
    
    leg2 = ax.legend(
        handles=patches2,
        loc='upper right',
        title="Municipalities of Charlevoix",
        fontsize=10,
        framealpha=1.0,
        bbox_to_anchor=(bbox1.x1, bbox1.y0 - LEGEND_GAP),
        bbox_transform=fig.transFigure,
        edgecolor='black',
        title_fontproperties={'weight': 'bold', 'size': 12}
    )
    ax.add_artist(leg2)
    
    fig.canvas.draw()
    bbox2 = leg2.get_window_extent().transformed(trans_inv)
    
    p0 = ax.transData.transform((0, 0))
    p1 = ax.transData.transform((10000, 0))
    fig_per_m = (
        (trans_inv.transform((p1[0] - p0[0], 0))[0] - trans_inv.transform((0, 0))[0])
        / 10000
    )
    
    ref_pops = [100, 1000, 10000]
    ref_labels = ['100', '1,000', '10,000']
    radii = [
        np.sqrt(p) * POPULATION_SCALE_FACTOR * fig_per_m
        for p in ref_pops
    ]
    
    leg3_h = (LEGEND_PAD + LEGEND_TITLE_HEIGHT + LEGEND_SPACING +
              max(radii) * 2 + LEGEND_SPACING + LEGEND_LABEL_HEIGHT + LEGEND_PAD)
    leg3_x = bbox2.x0
    leg3_y = bbox2.y0 - LEGEND_GAP - leg3_h
    leg3_w = bbox2.width
    
    ax.add_patch(mpatches.FancyBboxPatch(
        (leg3_x, leg3_y),
        leg3_w,
        leg3_h,
        boxstyle="round,pad=0,rounding_size=0.01",
        facecolor='white',
        edgecolor='black',
        linewidth=1,
        alpha=1.0,
        zorder=20,
        transform=fig.transFigure
    ))
    
    ax.text(
        leg3_x + leg3_w / 2,
        leg3_y + leg3_h - LEGEND_PAD,
        "Population Size",
        ha='center',
        va='top',
        fontsize=12,
        fontweight='bold',
        transform=fig.transFigure,
        zorder=21
    )
    
    col_w = (leg3_w - 2 * LEGEND_PAD) / 3
    circle_y = leg3_y + LEGEND_PAD + LEGEND_LABEL_HEIGHT + LEGEND_SPACING
    
    for i, (label, r) in enumerate(zip(ref_labels, radii)):
        cx = leg3_x + LEGEND_PAD + i * col_w + col_w / 2
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
            leg3_y + LEGEND_PAD + LEGEND_LABEL_HEIGHT / 2,
            label,
            ha='center',
            va='center',
            fontsize=10,
            transform=fig.transFigure,
            zorder=21
        )


# =============================================================================
# MAIN
# =============================================================================

def main():
    """Main execution."""
    print("=" * 60)
    print("SLSJ Genetic Contribution Map Generator")
    print("=" * 60)
    print()
    
    gdf = load_population_data("CharlevoixGradient.csv")
    genetic_df = load_genetic_contributions("CharlevoixGradientSubdivisions.csv")
    admin, lakes, rivers_ne, ocean, osm_rivers, osm_water = load_geodata_from_zip("geodata.zip")
    
    (gdf_proj, viewport_box, lakes_final, rivers_final, rivers_proj,
     ocean_final, osm_rivers_final, osm_water_final) = \
        project_and_clip_layers(gdf, lakes, rivers_ne, ocean, osm_rivers, osm_water)
    
    print("Aggregating by watercourse subdivision...")
    agg_data = aggregate_by_subdivision(gdf, gdf_proj, genetic_df)
    print(f"  ✓ Created {len(agg_data)} subdivision aggregations")
    
    # Calculate radii and apply dispersion
    print("Applying dispersion algorithm...")
    agg_data['radius_m'] = np.sqrt(agg_data['total_pop']) * POPULATION_SCALE_FACTOR
    agg_data['new_x'], agg_data['new_y'] = disperse_points(
        agg_data['centroid_x'].values,
        agg_data['centroid_y'].values,
        agg_data['radius_m'].values,
        agg_data
    )
    
    # Apply manual position adjustments
    print("Applying manual position adjustments...")
    for idx, row in agg_data.iterrows():
        subdivision = row['subdivision']
        if subdivision in MANUAL_POSITION_ADJUSTMENTS:
            x_offset, y_offset = MANUAL_POSITION_ADJUSTMENTS[subdivision]
            agg_data.at[idx, 'new_x'] += x_offset
            agg_data.at[idx, 'new_y'] += y_offset
            if x_offset != 0 or y_offset != 0:
                print(f"  {subdivision}: adjusted by ({x_offset:+.0f}, {y_offset:+.0f}) meters")
    
    print("Creating main map...")
    fig, ax = create_main_map(
        agg_data, lakes_final, rivers_final, ocean_final,
        osm_rivers_final, osm_water_final, rivers_proj, viewport_box,
        gdf_proj.total_bounds
    )
    
    print("Adding water body labels...")
    add_water_labels(ax, lakes_final, agg_data, osm_water_final)
    
    print("Adding legends...")
    add_legends(fig, ax, agg_data)
    
    print("Saving outputs...")
    plt.savefig("map_slsj_genetic.svg", dpi=300, format='svg', bbox_inches='tight')
    plt.savefig("map_slsj_genetic.png", dpi=300, format='png', bbox_inches='tight')
    
    print()
    print("=" * 60)
    print("✓ SUCCESS!")
    print("=" * 60)


if __name__ == "__main__":
    main()