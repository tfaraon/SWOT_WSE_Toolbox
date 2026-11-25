**fig01_synthetic**  
Controlled seven-pass test using the v0.3.2 synthetic flood at 100 m spacing. (A) Prescribed surface and estimates from a 1300 m window. The outlined square identifies the imposed peak removed if temporal rejection is enabled. (B) Number of contributing cells divided by the requested window area in cells. A missing elevation does not establish a dry surface.

**figS01_observatory_2025**  
Kati Thanda–Lake Eyre North. (A) Model bed elevations from the KTLE Observatory bathymetry cache, sampled on a 1.5 km display grid, and the two fixed stations in UTM zone 53S. The colour scale describes bed elevation in AHD. (B) Daily median WSE from the archived site export on 18 dates at Belt Bay and 21 at Madigan Gulf. Symbols are not interpolated. This archived series is distinct from the granule extraction in Fig. 2. Source: KTLE Observatory, revision df08a608a43afdcdfb0b927a071b40014c3dadc2, data/bathymetry.npz and site/data/swot_wse.json.

**fig02_station_levels**  
KTLE observations. Local surface elevation in the EGM2008 datum. Estimates are means of finite pixels in quality classes 0 and 1, in a nominal 1300 m square. Symbols distinguish grid spacing; each represents one accepted cycle/pass within a resolution. Dates are not interpolated. Windows at 100 and 250 m spacing span 1300 and 1250 m, respectively.

**fig03_sensitivity_100m**  
KTLE observations, 100 m grid. Left: retained/eligible cycle–pass counts for each configuration. Right: paired WSE differences from the pixel-screened 1300 m target window on passes accepted by every plotted method. Short black ticks indicate medians. The point, 500 m and 1300 m windows use per-cell quality screening. Centre QC gates the 1300 m window using the central cell. IQR rejection removes temporally flagged values from the 1300 m series.

**fig03_sensitivity_250m**  
KTLE observations, 250 m grid. Left: retained/eligible cycle–pass counts for each configuration. Right: paired WSE differences from the pixel-screened 1300 m target window on passes accepted by every plotted method. Short black ticks indicate medians. The point, 500 m and 1300 m windows use per-cell quality screening. Centre QC gates the 1300 m window using the central cell. IQR rejection removes temporally flagged values from the 1300 m series.

**fig04_footprint**  
KTLE observations. Belt Bay, 100 m raster. Illustrative dates target the first accepted date, maximum accepted WSE and maximum centroid offset. Repeated selections are filled by chronologically spaced dates. This selection is descriptive. Grid-cell centres are plotted relative to the fixed station. The diamond marks the centroid of cells contributing to each WSE estimate. The lower panels show WSE and two measures of spatial offset. Grey cells have unknown optical water status.