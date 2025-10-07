#!/usr/bin/env python
# coding: utf-8

import os
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
import warnings

warnings.filterwarnings('ignore')

def load_and_clean_data(ds, var_name):
    """Load and clean data by removing NaN values"""
    data = ds[var_name].values.flatten()
    return data[~np.isnan(data)]

def standardized_kde(data, x_grid):
    """Calculate KDE and normalize to have max value of 1"""
    kde = gaussian_kde(data)
    kde_values = kde(x_grid)
    return kde_values / kde_values.max()

def main():
    # Create output directory
    os.makedirs("../figs", exist_ok=True)
    
    print("Loading datasets...")
    ds_lgm = xr.open_dataset("../raw_data/lgmDA_lgm_Ocn_annualcropped.nc")
    ds_lh = xr.open_dataset("../raw_data/lgmDA_hol_Ocn_annualcropped.nc")
    
    # Load and clean data for all variables
    print("Processing data...")
    sst_lgm = load_and_clean_data(ds_lgm, 'sst')
    sst_lh = load_and_clean_data(ds_lh, 'sst')
    
    sss_lgm = load_and_clean_data(ds_lgm, 'sss')
    sss_lh = load_and_clean_data(ds_lh, 'sss')
    
    d18osw_lgm = load_and_clean_data(ds_lgm, 'd18osw')
    d18osw_lh = load_and_clean_data(ds_lh, 'd18osw')
    
    # Define colors for publication (colorblind-friendly)
    color_lgm = '#D55E00'  # Orange for LGM
    color_lh = '#0072B2'   # Blue for Late Holocene
    
    # Create figure with 2x3 subplots
    print("Creating figure...")
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    
    # Adjust spacing
    plt.subplots_adjust(hspace=0.35, wspace=0.3)
    
    # ============================================================================
    # Panel (a): SST KDE
    # ============================================================================
    ax = axes[0, 0]
    x_grid = np.linspace(min(sst_lgm.min(), sst_lh.min()), 
                         max(sst_lgm.max(), sst_lh.max()), 200)
    
    kde_lgm = standardized_kde(sst_lgm, x_grid)
    kde_lh = standardized_kde(sst_lh, x_grid)
    
    ax.plot(x_grid, kde_lgm, color=color_lgm, linewidth=2.5, label='LGM')
    ax.plot(x_grid, kde_lh, color=color_lh, linewidth=2.5, label='LH')
    ax.fill_between(x_grid, kde_lgm, alpha=0.3, color=color_lgm)
    ax.fill_between(x_grid, kde_lh, alpha=0.3, color=color_lh)
    
    ax.set_ylim(0, 1)
    ax.set_xlabel('SST [°C]', fontsize=13, fontweight='bold')
    ax.set_ylabel('Normalized Density', fontsize=13, fontweight='bold')
    ax.legend(loc='upper right', frameon=True, fontsize=11)
    ax.grid(True, alpha=0.3, linestyle='--', linewidth=0.5)
    ax.text(0.02, 0.98, '(a)', transform=ax.transAxes, 
            fontsize=14, fontweight='bold', va='top',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    ax.tick_params(labelsize=11)
    
    # ============================================================================
    # Panel (b): SSS KDE
    # ============================================================================
    ax = axes[0, 1]
    x_grid = np.linspace(min(sss_lgm.min(), sss_lh.min()), 
                         max(sss_lgm.max(), sss_lh.max()), 200)
    
    kde_lgm = standardized_kde(sss_lgm, x_grid)
    kde_lh = standardized_kde(sss_lh, x_grid)
    
    ax.plot(x_grid, kde_lgm, color=color_lgm, linewidth=2.5, label='LGM')
    ax.plot(x_grid, kde_lh, color=color_lh, linewidth=2.5, label='LH')
    ax.fill_between(x_grid, kde_lgm, alpha=0.3, color=color_lgm)
    ax.fill_between(x_grid, kde_lh, alpha=0.3, color=color_lh)
    
    ax.set_ylim(0, 1)
    ax.set_xlabel('SSS [PSU]', fontsize=13, fontweight='bold')
    ax.set_ylabel('Normalized Density', fontsize=13, fontweight='bold')
    ax.legend(loc='upper right', frameon=True, fontsize=11)
    ax.grid(True, alpha=0.3, linestyle='--', linewidth=0.5)
    ax.text(0.02, 0.98, '(b)', transform=ax.transAxes, 
            fontsize=14, fontweight='bold', va='top',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    ax.tick_params(labelsize=11)
    
    # ============================================================================
    # Panel (c): δ18Osw KDE
    # ============================================================================
    ax = axes[0, 2]
    x_grid = np.linspace(min(d18osw_lgm.min(), d18osw_lh.min()), 
                         max(d18osw_lgm.max(), d18osw_lh.max()), 200)
    
    kde_lgm = standardized_kde(d18osw_lgm, x_grid)
    kde_lh = standardized_kde(d18osw_lh, x_grid)
    
    ax.plot(x_grid, kde_lgm, color=color_lgm, linewidth=2.5, label='LGM')
    ax.plot(x_grid, kde_lh, color=color_lh, linewidth=2.5, label='LH')
    ax.fill_between(x_grid, kde_lgm, alpha=0.3, color=color_lgm)
    ax.fill_between(x_grid, kde_lh, alpha=0.3, color=color_lh)
    
    ax.set_ylim(0, 1)
    ax.set_xlabel('δ$^{18}$O$_{sw}$ [‰ VSMOW]', fontsize=13, fontweight='bold')
    ax.set_ylabel('Normalized Density', fontsize=13, fontweight='bold')
    ax.legend(loc='upper right', frameon=True, fontsize=11)
    ax.grid(True, alpha=0.3, linestyle='--', linewidth=0.5)
    ax.text(0.02, 0.98, '(c)', transform=ax.transAxes, 
            fontsize=14, fontweight='bold', va='top',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    ax.tick_params(labelsize=11)
    
    # ============================================================================
    # Panel (d): SST Boxplot
    # ============================================================================
    ax = axes[1, 0]
    
    bp = ax.boxplot([sst_lgm, sst_lh], 
                     labels=['LGM', 'LH'],
                     widths=0.6,
                     patch_artist=True,
                     showmeans=True,
                     meanprops=dict(marker='D', markerfacecolor='black', 
                                   markeredgecolor='black', markersize=6))
    
    # Color the boxes
    bp['boxes'][0].set_facecolor(color_lgm)
    bp['boxes'][0].set_alpha(0.7)
    bp['boxes'][1].set_facecolor(color_lh)
    bp['boxes'][1].set_alpha(0.7)
    
    # Style the box plot elements
    for element in ['whiskers', 'fliers', 'medians', 'caps']:
        plt.setp(bp[element], color='black', linewidth=1.5)
    
    ax.set_ylabel('SST [°C]', fontsize=13, fontweight='bold')
    ax.grid(True, alpha=0.3, linestyle='--', linewidth=0.5, axis='y')
    ax.text(0.02, 0.98, '(d)', transform=ax.transAxes, 
            fontsize=14, fontweight='bold', va='top',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    ax.tick_params(labelsize=11)
    
    # ============================================================================
    # Panel (e): SSS Boxplot
    # ============================================================================
    ax = axes[1, 1]
    
    bp = ax.boxplot([sss_lgm, sss_lh], 
                     labels=['LGM', 'LH'],
                     widths=0.6,
                     patch_artist=True,
                     showmeans=True,
                     meanprops=dict(marker='D', markerfacecolor='black', 
                                   markeredgecolor='black', markersize=6))
    
    # Color the boxes
    bp['boxes'][0].set_facecolor(color_lgm)
    bp['boxes'][0].set_alpha(0.7)
    bp['boxes'][1].set_facecolor(color_lh)
    bp['boxes'][1].set_alpha(0.7)
    
    # Style the box plot elements
    for element in ['whiskers', 'fliers', 'medians', 'caps']:
        plt.setp(bp[element], color='black', linewidth=1.5)
    
    ax.set_ylabel('SSS [PSU]', fontsize=13, fontweight='bold')
    ax.grid(True, alpha=0.3, linestyle='--', linewidth=0.5, axis='y')
    ax.text(0.02, 0.98, '(e)', transform=ax.transAxes, 
            fontsize=14, fontweight='bold', va='top',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    ax.tick_params(labelsize=11)
    
    # ============================================================================
    # Panel (f): δ18Osw Boxplot
    # ============================================================================
    ax = axes[1, 2]
    
    bp = ax.boxplot([d18osw_lgm, d18osw_lh], 
                     labels=['LGM', 'LH'],
                     widths=0.6,
                     patch_artist=True,
                     showmeans=True,
                     meanprops=dict(marker='D', markerfacecolor='black', 
                                   markeredgecolor='black', markersize=6))
    
    # Color the boxes
    bp['boxes'][0].set_facecolor(color_lgm)
    bp['boxes'][0].set_alpha(0.7)
    bp['boxes'][1].set_facecolor(color_lh)
    bp['boxes'][1].set_alpha(0.7)
    
    # Style the box plot elements
    for element in ['whiskers', 'fliers', 'medians', 'caps']:
        plt.setp(bp[element], color='black', linewidth=1.5)
    
    ax.set_ylabel('δ$^{18}$O$_{sw}$ [‰ VSMOW]', fontsize=13, fontweight='bold')
    ax.grid(True, alpha=0.3, linestyle='--', linewidth=0.5, axis='y')
    ax.text(0.02, 0.98, '(f)', transform=ax.transAxes, 
            fontsize=14, fontweight='bold', va='top',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    ax.tick_params(labelsize=11)
    
    # Save figure
    print("Saving figure...")
    plt.savefig("../figs/kde_boxplot_comparison.png", dpi=300, bbox_inches="tight")
    plt.close()
    
    print("\n" + "=" * 80)
    print("PROCESSING COMPLETE")
    print("=" * 80)
    print("\nFigure saved to ../figs/")
    print("  - kde_boxplot_comparison.png")
    print("=" * 80)


if __name__ == "__main__":
    main()
