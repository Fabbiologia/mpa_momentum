#!/usr/bin/env python3
"""Rasterise UNEP-WCMC marine habitat layers at 1 km global resolution.

This script expects the raw shapefiles (points and/or polygons) to live in
``data/raw/habitats_original``.  For each habitat we:
    * read all shapefiles that belong to the habitat (e.g. points + polygons)
    * project geometries to the equal-area EPSG:6933 CRS
    * buffer point data to 1 km² circles so that all geometries represent area
    * clean invalid geometries and dissolve overlaps
    * rasterise to a 1 km grid covering the habitat extent
    * write a GeoTIFF to ``data/processed/habitat_rasters/<habitat>_1km.tif``
    * compute total area (km²) covered by the raster and store a summary CSV

The output rasters use pixel values of 1 for habitat presence (0 for absence).
The pipeline follows Kumagai et al. (2021) in using mapped UNEP-WCMC layers
instead of modelled surfaces.
"""

from __future__ import annotations

import math
import os
import re
from pathlib import Path
from typing import Dict, List

import geopandas as gpd
import numpy as np
import pandas as pd
import rasterio
from rasterio import features
from rasterio.transform import from_origin
from shapely.geometry import Point
from shapely.ops import unary_union

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
PROJECT_ROOT = Path(__file__).resolve().parents[1]
RAW_HABITAT_DIR = PROJECT_ROOT / "data" / "raw" / "habitats_original"
OUTPUT_RASTER_DIR = PROJECT_ROOT / "data" / "processed" / "habitat_rasters"
OUTPUT_RASTER_DIR.mkdir(parents=True, exist_ok=True)
SUMMARY_PATH = PROJECT_ROOT / "data" / "processed" / "habitat_raster_summary.csv"

TARGET_CRS = "EPSG:6933"  # Equal-area CRS (Equal Earth)
RESOLUTION = 1000  # metres (≈ 1 km)
BUFFER_RADIUS_METRES = math.sqrt(1 / math.pi) * 1000  # 1 km² circle radius

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _clean_habitat_name(filename: str) -> str:
    """Infer a habitat key from the UNEP-WCMC file names."""
    name = Path(filename).stem.lower()
    name = re.sub(r"^wcmc\d+_?", "", name)
    name = re.sub(r"_v\d+_\d+$", "", name)

    if "mangrove" in name or "gmw" in name:
        return "mangroves"
    if "seagrass" in name:
        return "seagrasses"
    if "saltmarsh" in name:
        return "saltmarshes"
    if "coralreef" in name or "coral_reef" in name:
        return "coral_reefs"
    if "coldcoral" in name or "cold_corals" in name:
        return "cold_water_corals"
    if "knolls" in name or "seamount" in name:
        return "seamounts_knolls"
    return name


def _group_habitat_files() -> Dict[str, List[Path]]:
    if not RAW_HABITAT_DIR.exists():
        raise FileNotFoundError(
            f"Habitat directory {RAW_HABITAT_DIR} does not exist. "
            "Please place the UNEP-WCMC shapefiles there."
        )

    groups: Dict[str, List[Path]] = {}
    for shp in RAW_HABITAT_DIR.glob("*.shp"):
        if "metadata" in shp.name.lower():
            continue
        key = _clean_habitat_name(shp.name)
        groups.setdefault(key, []).append(shp)
    return groups


def _read_and_prepare(path: Path) -> gpd.GeoDataFrame:
    gdf = gpd.read_file(path)
    if gdf.empty:
        return gdf
    if gdf.crs is None:
        raise ValueError(f"Dataset {path} is missing a CRS definition")
    if gdf.crs.to_string() != TARGET_CRS:
        gdf = gdf.to_crs(TARGET_CRS)
    return gdf


def _buffer_points(gdf: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    buffered = gdf.copy()
    buffered["geometry"] = buffered.geometry.apply(
        lambda geom: geom.buffer(BUFFER_RADIUS_METRES) if isinstance(geom, Point) else geom
    )
    return buffered


def _dissolve_geometries(geoms: gpd.GeoSeries) -> gpd.GeoSeries:
    union = unary_union(geoms)
    dissolved = gpd.GeoSeries([union], crs=TARGET_CRS)
    dissolved[0] = dissolved.iloc[0].buffer(0)
    return dissolved


def _rasterise(habitat: str, geom_series: gpd.GeoSeries) -> Dict[str, float]:
    bounds = geom_series.total_bounds
    minx, miny, maxx, maxy = bounds
    width = int(np.ceil((maxx - minx) / RESOLUTION))
    height = int(np.ceil((maxy - miny) / RESOLUTION))

    transform = from_origin(minx, maxy, RESOLUTION, RESOLUTION)

    shapes = [(geom, 1) for geom in geom_series.geometry if not geom.is_empty]

    raster = features.rasterize(
        shapes,
        out_shape=(height, width),
        transform=transform,
        fill=0,
        dtype="uint8"
    )

    raster_path = OUTPUT_RASTER_DIR / f"{habitat}_1km.tif"
    with rasterio.open(
        raster_path,
        "w",
        driver="GTiff",
        height=height,
        width=width,
        count=1,
        dtype="uint8",
        crs=TARGET_CRS,
        transform=transform,
        compress="lzw",
        tiled=True,
        blockxsize=256,
        blockysize=256,
    ) as dst:
        dst.write(raster, 1)

    cell_area_km2 = (RESOLUTION * RESOLUTION) / 1e6
    area_km2 = raster.sum() * cell_area_km2

    return {"habitat": habitat, "cells": int(raster.sum()), "area_km2": float(area_km2)}


def main():
    groups = _group_habitat_files()
    if not groups:
        raise RuntimeError(
            "No habitat shapefiles found. Please ensure data/raw/habitats_original "
            "contains the UNEP-WCMC layers (.shp files)."
        )

    summary_rows: List[Dict[str, float]] = []

    for habitat, paths in sorted(groups.items()):
        print(f"Processing {habitat} ({len(paths)} files)")

        frames = []
        for path in paths:
            gdf = _read_and_prepare(path)
            if gdf.empty:
                continue
            if gdf.geom_type.isin(["Point", "MultiPoint"]).all():
                print("  buffering point features")
                gdf = _buffer_points(gdf)
            frames.append(gdf)

        if not frames:
            print("  warning: no geometries found; skipping")
            continue

        merged = gpd.GeoDataFrame(pd.concat(frames, ignore_index=True), crs=TARGET_CRS)
        merged["geometry"] = merged.geometry.buffer(0)

        dissolved = _dissolve_geometries(merged.geometry)
        stats = _rasterise(habitat, dissolved)
        summary_rows.append(stats)

    summary = pd.DataFrame(summary_rows).sort_values("habitat")
    summary.to_csv(SUMMARY_PATH, index=False)
    print(f"Summary written to {SUMMARY_PATH}")


if __name__ == "__main__":
    main()
