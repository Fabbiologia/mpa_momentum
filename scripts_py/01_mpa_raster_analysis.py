#!/usr/bin/env python3
"""
Rasterise ProtectedSeas MPAs at 1 km resolution and compute annual LFP totals.

Outputs:
    - <output-dir>/lfp_max_1km.tif : single-band GeoTIFF with max LFP value (0-5).
    - <output-dir>/lfp_yearly_counts.csv : yearly counts of 1 km cells per LFP value.
"""

import argparse
import multiprocessing as mp
import time
import warnings
from datetime import datetime
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import rasterio
from joblib import Parallel, delayed
from rasterio.enums import MergeAlg
from rasterio.features import rasterize
from rasterio.transform import from_bounds
from tqdm import tqdm

warnings.filterwarnings("ignore")


def load_datasets(data_dir: Path):
    shp_path = data_dir / "Navigator_AllSites_032825_shp" / "Navigator_AllSites_032825.shp"
    csv_candidates = [
        data_dir / "Navigator_AllSites_032825_corrected.csv",
        data_dir / "Navigator_AllSites_032825.csv",
    ]
    csv_path = next((p for p in csv_candidates if p.exists()), None)
    if csv_path is None:
        raise FileNotFoundError("No Navigator_AllSites CSV file found in provided data directory.")

    print(f"Reading shapefile from {shp_path}")
    gdf = gpd.read_file(shp_path)
    print(f"Shapefile loaded with {len(gdf)} records")

    print(f"Reading metadata from {csv_path}")
    metadata = pd.read_csv(csv_path, low_memory=False)
    print(f"Metadata loaded with {len(metadata)} records")
    return gdf, metadata


def join_attributes(gdf: gpd.GeoDataFrame, metadata: pd.DataFrame) -> gpd.GeoDataFrame:
    if "SITE_ID" in gdf.columns and "site_id" in metadata.columns:
        gdf = gdf.copy()
        gdf["site_id"] = gdf["SITE_ID"]
        merged = gdf.merge(metadata, on="site_id", how="left")
    else:
        common_fields = set(gdf.columns).intersection(metadata.columns)
        join_field = next(
            (field for field in common_fields if any(key in field.lower() for key in ("id", "code", "key"))),
            None,
        )
        if join_field is None and common_fields:
            join_field = list(common_fields)[0]
        if join_field is None:
            print("Warning: could not identify join field; returning original geometries.")
            return gdf
        merged = gdf.merge(metadata, on=join_field, how="left")
    print(f"Merged dataset has {len(merged)} records")
    return merged


def prepare_subset(
    merged: gpd.GeoDataFrame,
    category_mode: str = "mpa",
    exclude_categories=None,
    include_high_seas: bool = True,
) -> gpd.GeoDataFrame:
    subset = merged.copy()

    if "category_name" in subset.columns:
        if category_mode == "mpa":
            subset = subset[subset["category_name"] == "Marine Protected Area"].copy()
        elif category_mode == "other":
            subset = subset[subset["category_name"] != "Marine Protected Area"].copy()
        print(
            f"Subset after category filter ({category_mode}): {len(subset)} records"
        )
    else:
        print("Warning: 'category_name' column missing; using all geometries.")

    if exclude_categories:
        before = len(subset)
        subset = subset[~subset["category_name"].isin(exclude_categories)].copy()
        removed = before - len(subset)
        if removed > 0:
            print(f"Excluded {removed} records via category exclusion list.")

    if not include_high_seas and "country" in subset.columns:
        before = len(subset)
        subset = subset[subset["country"] != "High Seas / International"].copy()
        removed = before - len(subset)
        if removed > 0:
            print(f"Removed {removed} high seas / international records.")

    subset["lfp"] = pd.to_numeric(subset["lfp"], errors="coerce").fillna(0)

    if hasattr(subset.geometry, "make_valid"):
        subset = subset.set_geometry(subset.geometry.make_valid())
    else:
        subset = subset.set_geometry(subset.buffer(0))

    year_numeric = pd.to_numeric(subset["year_est"], errors="coerce")
    year_valid = year_numeric.where((year_numeric >= 1800) & (year_numeric <= 2025))
    subset["year_est"] = year_valid.fillna(1999).astype(int)
    return subset


def rasterize_lfp(mpas: gpd.GeoDataFrame, resolution: int, output_raster: Path):
    target_crs = "EPSG:6933"
    mpas_projected = mpas.to_crs(target_crs)

    bounds = mpas_projected.total_bounds
    width = int(np.ceil((bounds[2] - bounds[0]) / resolution))
    height = int(np.ceil((bounds[3] - bounds[1]) / resolution))
    transform = from_bounds(bounds[0], bounds[1], bounds[2], bounds[3], width, height)

    print(f"Raster grid: {width:,} x {height:,} ({resolution} m resolution)")

    shapes = []
    for _, row in tqdm(
        mpas_projected.iterrows(),
        total=len(mpas_projected),
        desc="Preparing shapes",
        unit="feature",
    ):
        if row.geometry is not None and not row.geometry.is_empty:
            try:
                lfp_value = float(row.lfp)
                shapes.append((row.geometry, lfp_value))
            except (ValueError, TypeError):
                continue

    lfp_raster = np.zeros((height, width), dtype=np.float32)
    unique_lfp = sorted({value for _, value in shapes}, reverse=True)
    print(f"Rasterising for LFP values: {unique_lfp}")

    for lfp_value in tqdm(unique_lfp, desc="Rasterising by LFP"):
        shapes_with_lfp = [(geom, 1) for geom, value in shapes if value == lfp_value]
        if not shapes_with_lfp:
            continue

        mask = rasterize(
            shapes_with_lfp,
            out_shape=(height, width),
            transform=transform,
            fill=0,
            all_touched=True,
            dtype=np.uint8,
            merge_alg=MergeAlg.add,
        )

        lfp_raster = np.where((mask > 0) & (lfp_raster == 0), lfp_value, lfp_raster)

    output_raster.parent.mkdir(parents=True, exist_ok=True)
    with rasterio.open(
        output_raster,
        "w",
        driver="GTiff",
        height=height,
        width=width,
        count=1,
        dtype=np.float32,
        crs=target_crs,
        transform=transform,
        compress="lzw",
    ) as dst:
        dst.write(lfp_raster.astype(np.float32), 1)

    print(f"Saved LFP raster to {output_raster}")
    return mpas_projected, transform, (height, width)


def annual_lfp_totals(
    mpas_projected: gpd.GeoDataFrame,
    transform,
    shape,
    resolution: int,
    start_year: int,
    end_year: int,
    n_jobs: int,
):
    height, width = shape
    years = list(range(start_year, end_year + 1))

    def process_year(year_idx, year):
        subset = mpas_projected[mpas_projected["year_est"] <= year]
        year_totals = {f"lfp_{i}": 0 for i in range(1, 6)}
        for lfp_value in range(1, 6):
            geometries = subset[subset["lfp"] == lfp_value].geometry
            shapes = [(geom, 1) for geom in geometries if geom is not None and not geom.is_empty]
            if not shapes:
                continue
            mask = rasterize(
                shapes,
                out_shape=(height, width),
                transform=transform,
                fill=0,
                all_touched=True,
                dtype=np.uint8,
                merge_alg=MergeAlg.add,
            )
            year_totals[f"lfp_{lfp_value}"] = int(mask.astype(bool).sum())
        return year_idx, year_totals

    results = Parallel(n_jobs=n_jobs)(
        delayed(process_year)(idx, year) for idx, year in enumerate(years)
    )

    totals = {"year": years}
    for lfp_value in range(1, 6):
        totals[f"lfp_{lfp_value}"] = [0] * len(years)

    for year_idx, year_totals in results:
        for key, value in year_totals.items():
            totals[key][year_idx] = value

    return pd.DataFrame(totals)


def main():
    parser = argparse.ArgumentParser(description="Rasterise MPAs and compute annual LFP totals.")
    parser.add_argument("--data-dir", type=Path, default=Path("data/raw/protected_seas_navigator_all_sites"))
    parser.add_argument("--output-dir", type=Path, default=Path("data/processed/mpa_raster"))
    parser.add_argument("--resolution", type=int, default=1000, help="Raster resolution in metres.")
    parser.add_argument("--start-year", type=int, default=2000)
    parser.add_argument("--end-year", type=int, default=datetime.now().year)
    parser.add_argument("--max-cores", type=int, default=10)
    parser.add_argument(
        "--category-mode",
        choices=["mpa", "other", "all"],
        default="mpa",
        help="Filter geometries by category (MPAs only, non-MPAs, or all).",
    )
    parser.add_argument(
        "--exclude-categories",
        nargs="*",
        default=None,
        help="category_name values to exclude after initial filtering.",
    )
    parser.add_argument(
        "--include-high-seas",
        action="store_true",
        help="Keep records tagged as 'High Seas / International'.",
    )
    parser.add_argument(
        "--output-prefix",
        type=str,
        default=None,
        help="Prefix for output filenames (defaults to the category mode).",
    )
    args = parser.parse_args()

    start_time = time.time()
    data_dir = args.data_dir
    output_dir = args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    gdf, metadata = load_datasets(data_dir)
    merged = join_attributes(gdf, metadata)

    exclude_categories = args.exclude_categories or []
    if args.category_mode == "other" and not exclude_categories:
        exclude_categories = ["Fisheries Management Area"]

    subset = prepare_subset(
        merged,
        category_mode=args.category_mode,
        exclude_categories=exclude_categories,
        include_high_seas=args.include_high_seas,
    )

    if subset.empty:
        raise RuntimeError("No geometries remain after applying filters.")

    num_cores = mp.cpu_count()
    n_jobs = min(args.max_cores, num_cores)
    print(f"Using {n_jobs} worker processes")

    prefix = args.output_prefix or args.category_mode
    raster_path = output_dir / f"{prefix}_lfp_max_1km.tif"
    subset_projected, transform, shape = rasterize_lfp(subset, args.resolution, raster_path)

    print("Computing annual LFP totals (cell counts)...")
    totals_df = annual_lfp_totals(
        subset_projected,
        transform,
        shape,
        args.resolution,
        args.start_year,
        args.end_year,
        n_jobs,
    )
    totals_path = output_dir / f"{prefix}_lfp_yearly_counts.csv"
    totals_df.to_csv(totals_path, index=False)
    print(f"Saved annual totals to {totals_path}")

    elapsed = time.time() - start_time
    print(f"Processing completed in {elapsed:.1f} seconds")


if __name__ == "__main__":
    main()
