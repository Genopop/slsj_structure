import csv
import json
import pandas as pd
import pickle
import geneakit as gen

# Get paths to data files
with open("../paths.json", 'r') as file:
    paths = json.load(file)

# Load the BALSAC genealogy
ped = gen.genealogy(paths['balsac_genealogy'])

# Create dictionaries for municipality and region codes
city_code_to_string = {}
with open(paths['geography_definitions'], 'r') as file:
    reader = csv.reader(file)
    next(reader)  # Skip header row
    for row in reader:
        UrbMariage, UrbIdMariage, *_ = row
        city_code_to_string[int(UrbIdMariage)] = UrbMariage
region_code_to_string = {
    25705: "Abitibi",
    25706: "Bas-Saint-Laurent",
    25707: "Beauce",
    25708: "Bois-Francs",
    25709: "Charlevoix",
    25710: "Côte-de-Beaupré",
    25711: "Côte-du-Sud",
    25712: "Côte-Nord",
    25713: "Estrie",
    25714: "Gaspésie",
    25715: "Île-de-Montréal",
    25716: "Îles-de-la-Madeleine",
    25717: "Lanaudière",
    25718: "Laurentides",
    25719: "Mauricie",
    25720: "Outaouais",
    25721: "Agglomération de Québec",
    25722: "Région de Québec",
    25723: "Nord du Québec",
    25724: "Richelieu",
    25725: "Rive Nord-Ouest de Montréal",
    25726: "Rive Sud de Montréal",
    25727: "Saguenay–Lac-Saint-Jean",
    25728: "Témiscamingue",
    27118: "Côte-de-Beaupré",
    27119: "Portneuf",
    27120: "Lévis-Lotbinière"
}

# Create dictionaries for proband and parent information
city_proband = {}
region_proband = {}
year_proband = {}
city_parent = {}
region_parent = {}
year_parent = {}

with open(paths['demography_information'], 'r', encoding='cp1252') as file:
    reader = csv.reader(file, delimiter=';')
    next(reader)  # Skip header row
    for row in reader:
        if len(row) == 0: break
        IndID, CaG, ERRQ, PereID, MereID, Sexe, PaysOrigine, DateNaissance, RegionNaissance, \
        DateDeces, RegionDeces, DateMariage, ConjointID, URBMariage, RegionMariage, \
        DateMariageParents, URBMariageParents, RegionMariageParents = row

        IndID = int(IndID) if IndID != 'NA' else 0
        URBMariage = int(URBMariage) if URBMariage != 'NA' else 0
        RegionMariage = int(RegionMariage) if RegionMariage != 'NA' else 0
        DateMariage = int(DateMariage) if DateMariage != 'NA' else 0
        URBMariageParents = int(URBMariageParents) if URBMariageParents != 'NA' else 0
        RegionMariageParents = int(RegionMariageParents) if RegionMariageParents != 'NA' else 0
        DateMariageParents = int(DateMariageParents) if DateMariageParents != 'NA' else 0

        city_proband[IndID] = city_code_to_string.get(URBMariage, 0)
        region_proband[IndID] = region_code_to_string.get(RegionMariage, 0)
        year_proband[IndID] = DateMariage
        city_parent[IndID] = city_code_to_string.get(URBMariageParents, 0)
        region_parent[IndID] = region_code_to_string.get(RegionMariageParents, 0)
        year_parent[IndID] = DateMariageParents

# Extract the probands who married in Saguenay–Lac-Saint-Jean between 1931 and 1960
inds_all = sorted(ped.keys())
inds_in_slsj = [ind for ind in inds_all if region_proband.get(ind, 'Inconnue') == 'Saguenay–Lac-Saint-Jean']
inds_from_1931 = [ind for ind in inds_all if year_proband.get(ind, 0) >= 1935]
inds_until_1960 = [ind for ind in inds_all if year_proband.get(ind, 1961) <= 1960]
set_slsj = set(inds_in_slsj)
set_1931 = set(inds_from_1931)
set_1960 = set(inds_until_1960)
inds_intersection = [ind for ind in inds_all if ind in set_slsj and ind in set_1931 and ind in set_1960]

# Create a new genealogy that starts from those probands
iso_ped = gen.branching(ped, pro=inds_intersection)

# Remove parents from the probands
pro = gen.pro(iso_ped)

# Define the watercourse subdivisions
city_to_water_boundary = {
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
    "Saint-Prime": "Ashuapmushuan-Mistassini",
    "Saint-Stanislas-du-Lac-Saint-Jean": "Mistassini-Péribonka",
    "Saint-Thomas-Didyme": "Ashuapmushuan-Mistassini",
    "Sainte-Hedwidge": "Métabetchouane-Ashuapmushuan",
    "Sainte-Jeanne-d'Arc-du-Lac-Saint-Jean": "Mistassini-Péribonka",
    "Sainte-Monique-de-Honfleur": "East of Lac Saint-Jean",
    "Sainte-Rose-du-Nord": "North of Saguenay River",
    "Sainte-Élisabeth-de-Proulx": "Mistassini-Péribonka",
    "Shipshaw": "North of Saguenay River",
    "Val-Jalbert": "Métabetchouane-Ashuapmushuan",
    'Baie-Sainte-Catherine': 'Charlevoix',
    'Baie-Saint-Paul': 'Charlevoix',
    "Cap A L'aigle": 'Charlevoix',
    'Clermont': 'Charlevoix',
    'La Malbaie': 'Charlevoix',
    'Les Éboulements': 'Charlevoix',
    'Notre Dame Des Monts': 'Charlevoix',
    'Petite-Rivière-Saint-François': 'Charlevoix',
    'Pointe-au-Pic': 'Charlevoix',
    'St Aime Des Lacs': 'Charlevoix',
    "Saint-Bernard-de-l'Isle-aux-Coudres": 'Charlevoix',
    'Sainte-Agnès': 'Charlevoix',
    'Saint-Fidèle-de-Mont-Murray': 'Charlevoix',
    'Saint-Hilarion': 'Charlevoix',
    'Saint-Irénée': 'Charlevoix',
    'Saint-Joseph-de-la-Rive': 'Charlevoix',
    "Saint-Louis-de-l'Isle-aux-Coudres": 'Charlevoix',
    'Saint-Siméon': 'Charlevoix',
    'Saint-Urbain-de-Charlevoix': 'Charlevoix'
}

# Initiate the DataFrame for the founders' genetic contributions
founder_contributions_df = pd.DataFrame(
    columns=['Founder ID', 'Year', 'Region', 'Origin',
             'Individual GC', 'Watercourse', 'Municipality'])

# Recursively add genetic contributions to the DataFrame
def add_founder_contribution(individual, year, depth, proband):
    if region_parent.get(individual.ind, 'Unknown') != 'Saguenay–Lac-Saint-Jean':
        founder = individual.ind
        region = region_parent[individual.ind]
        origin = city_parent[individual.ind]
        contribution = 0.5 ** depth
        watercourse = city_to_water_boundary[city_proband[proband]]
        municipality = city_proband[proband]
        founder_contributions_df.loc[founder_contributions_df.shape[0]] = [
            founder, year, region, origin,
            contribution, watercourse, municipality]
    else:
        if individual.father.ind != 0:
            add_founder_contribution(individual.father, year_parent[individual.ind], depth+1, proband)
        if individual.mother.ind != 0:
            add_founder_contribution(individual.mother, year_parent[individual.ind], depth+1, proband)

# Start from the probands and move up
for ind in pro:
	add_founder_contribution(ped[ind], year_proband[ind], 0, ind)
	
with open(paths['wd'] + "results/pickles/balsac_contributions.pkl", 'wb') as file:
    pickle.dump(founder_contributions_df, file)