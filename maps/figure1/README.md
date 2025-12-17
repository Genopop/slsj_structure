# Saguenay–Lac-Saint-Jean Population Distribution Map

## Overview

This repository contains the code to generate a cartographic visualization of population distribution across the Saguenay–Lac-Saint-Jean region of Quebec, Canada. The map displays municipalities as proportionally-sized circles colored by watercourse subdivision, with geographic features from Natural Earth and OpenStreetMap.

## Requirements

### Python Version
- Python 3.8 or higher

### Dependencies
Install required packages using:
```bash
pip install -r requirements.txt
```

Key dependencies:
- `geopandas` - Geospatial data handling
- `matplotlib` - Visualization
- `shapely` - Geometric operations
- `numpy`, `pandas`, `scipy` - Numerical computing
- `requests` - API queries

## Input Data

The script requires a CSV file named `CharlevoixGradient.csv` with the following columns:

- `Longitude` - Decimal degrees (WGS84)
- `Latitude` - Decimal degrees (WGS84)
- `Location of proband's marriage` - Municipality name
- `Population size` - Population count (optional, defaults to 1000)

## Usage

### Basic Usage
```bash
python map_slsj.py
```

This will:
1. Load population data from `CharlevoixGradient.csv`
2. Download geographic data from Natural Earth
3. Query OpenStreetMap for detailed river information
4. Generate the map
5. Save outputs as `map_slsj.svg` and `map_slsj.png`

### Expected Runtime
- First run: 2-5 minutes (downloading geographic data)
- Subsequent runs: 1-2 minutes (if geographic data cached)

## Output Files

- `map_slsj.svg` - Vector format (recommended for publications)
- `map_slsj.png` - Raster format at 300 DPI

## Configuration

All visualization parameters are defined as constants at the top of the script:

### Key Parameters

```python
# Population circle sizing
POPULATION_SCALE_FACTOR = 40.0      # Adjust circle sizes
LABEL_SIZE_THRESHOLD = 2500          # Minimum population for city labels

# Spatial extent
BUFFER_METERS = 25000                # Map margin in meters

# Circle positioning
DISPERSION_ITERATIONS = 50           # Overlap adjustment iterations
REPULSION_FACTOR = 0.5               # Strength of circle repulsion
ATTRACTION_FACTOR = 0.05             # Attraction to original position
```

### Color Scheme

Watercourse subdivisions are defined in `SUBDIVISION_COLORS`:
```python
SUBDIVISION_COLORS = {
    "East of Ha! Ha!": 'gold',
    "South of Saguenay River": 'lime',
    "North of Saguenay River": 'deeppink',
    # ... etc
}
```

## Algorithm Details

### Circle Dispersion
The script uses a force-directed layout algorithm to prevent circle overlap while minimizing displacement from original geographic positions. This runs for 50 iterations by default.

### River Labeling
Text follows river paths using spline interpolation with:
- Geometry smoothing to reduce noise
- Angular smoothing for readable character orientation
- Configurable placement along specific river segments

### Coordinate Systems
- Input data: WGS84 (EPSG:4326)
- Map projection: NAD83 / Quebec Lambert (EPSG:32198)

## Data Sources

1. **Natural Earth** (public domain)
   - Administrative boundaries (1:50m)
   - Lakes (1:10m)
   - Rivers (1:10m)
   - Ocean (1:10m)

2. **OpenStreetMap** (ODbL license)
   - Detailed river geometries
   - Water body polygons

## Troubleshooting

### Common Issues

**Import errors**
```bash
pip install --upgrade geopandas matplotlib shapely
```

**OpenStreetMap query timeout**
- Increase `OSM_TIMEOUT` constant (default: 180 seconds)
- Check internet connection

**Missing data file**
```
FileNotFoundError: CharlevoixGradient.csv
```
- Ensure CSV file is in working directory
- Check column names match expected format

**Memory errors with large datasets**
- Reduce `BUFFER_METERS` to limit geographic extent
- Decrease `DISPERSION_ITERATIONS`

## Author

- Gilles-Philippe Morin, with the help of Gemini 3 Pro and Claude Sonnet 4.5