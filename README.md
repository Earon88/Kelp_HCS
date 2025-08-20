Spatio-temporal dynamics of giant kelp forests in the Humboldt Current System

Dataset DOI: 10.5061/dryad.vt4b8gv58

Description of the data and file structure

The database represents the extent of the floating canopy of M. pyrifera along the Chilean part of the Humboldt Current System and its seasonal variation from summer 2016 to summer 2015.

Variable Name: Data Type;
Unit/Format
Description

season: Categorical (Factor); 
Text (Fall, Winter, Spring, Summer) 
The season of the year when the observation was recorded. 

year: Integer;
YYYY
The year when the observation was recorded.

area: Integer 
meters squares 
A unique code identifying a specific maritime area or zone.


ECOREGION: Categorical (Factor);
Text (Araucanian, Central Chile, Humboldtian)
The marine ecoregion where the observation was recorded, based on the Marine Ecoregions of the World (MEOW) framework.


quarter Character
Text ("Q1", "Q2", "Q3", "Q4")
The quarter of the year (Q1, Q2, Q3, Q4) when the observation was recorded.

yearqtr: Numeric (Yearqtr)
YYYY.q
A combined representation of the year and quarter for time series analysis.


lat: Numeric
Decimal degrees (WGS84)
The latitude coordinate

lon: Numeric
Decimal degrees (WGS84)
The longitude coordinate

X: Numeric
Easting Kilometers EPSG:32719
The Easting coordinate in a projected coordinate system (likely UTM Zone 19S, EPSG:32719).

Y: Numeric
Northing Kilometers EPSG:32719
The Northing coordinate in a projected coordinate system (likely UTM Zone 19S, EPSG:32719).


Files and variables

File: Kelp_detection_database.R

Description: Code for generating polygons of marine ecoregions and processing data obtained from image processing in Google Earth Engine

File: MS_ST_floating_canopy.R

Description: Code to analyze the spatiotemporal database of the floating canopy of M. pyrifera

File: MP_data

Description: Data in format RDS, used to model the spatiotemporal variation of the floating canopy of M. pyrifera

File: MP_data.csv

Description: Data in format csv, used to model the spatiotemporal variation of the floating canopy of M. pyrifera

Code/software

Data collection and processing was carried out using R and Google Earth Engine.

https://code.earthengine.google.com/fc06f412d7b0d57b53765e29d82e11e2?hl=es-419


Data was derived from the following sources:
Sentinel-2
