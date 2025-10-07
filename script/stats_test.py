#!/usr/bin/env python
"""
Rigorous Statistical Analysis: LGM vs Late Holocene
Comprehensive non-parametric testing with multiple normality assessments,
effect sizes, and complexity measures.
"""

import pandas as pd
import numpy as np
from scipy import stats
from scipy.stats import entropy
import os
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')


def cliffs_delta(x, y):
    """
    Calculate Cliff's Delta effect size.
    Returns: delta value and interpretation
    """
    nx, ny = len(x), len(y)
    if nx == 0 or ny == 0:
        return np.nan, "undefined"
    
    # Count dominances
    dominance = sum(1 if xi > yi else (-1 if xi < yi else 0) 
                    for xi in x for yi in y)
    delta = dominance / (nx * ny)
    
    # Interpretation
    abs_delta = abs(delta)
    if abs_delta < 0.147:
        interp = "negligible"
    elif abs_delta < 0.330:
        interp = "small"
    elif abs_delta < 0.474:
        interp = "medium"
    else:
        interp = "large"
    
    return delta, interp


def approximate_entropy(x, m=2, r=None):
    """Calculate Approximate Entropy for complexity measure"""
    x = np.array(x)
    n = len(x)
    
    if r is None:
        r = 0.2 * np.std(x)
    
    def _maxdist(xi, xj):
        return max([abs(ua - va) for ua, va in zip(xi, xj)])
    
    def _phi(m):
        patterns = np.array([[x[j] for j in range(i, i + m)] 
                            for i in range(n - m + 1)])
        C = np.zeros(n - m + 1)
        for i in range(n - m + 1):
            C[i] = sum(1 for j in range(n - m + 1) 
                      if _maxdist(patterns[i], patterns[j]) <= r) / (n - m + 1)
        return np.sum(np.log(C)) / (n - m + 1)
    
    return abs(_phi(m) - _phi(m + 1))


def sample_entropy(x, m=2, r=None):
    """Calculate Sample Entropy for complexity measure"""
    x = np.array(x)
    n = len(x)
    
    if r is None:
        r = 0.2 * np.std(x)
    
    def _maxdist(xi, xj):
        return max([abs(ua - va) for ua, va in zip(xi, xj)])
    
    def _phi(m):
        patterns = np.array([[x[j] for j in range(i, i + m)] 
                            for i in range(n - m + 1)])
        count = 0
        for i in range(n - m):
            for j in range(i + 1, n - m + 1):
                if _maxdist(patterns[i], patterns[j]) <= r:
                    count += 1
        return count
    
    A = _phi(m + 1)
    B = _phi(m)
    
    if B == 0:
        return np.nan
    
    return -np.log(A / B)


def perform_normality_tests(data, var_name, group_name):
    """Perform multiple normality tests"""
    results = {}
    
    # Remove NaN values
    clean_data = data[~np.isnan(data)]
    n = len(clean_data)
    
    if n < 3:
        return {"error": f"Insufficient data (n={n})"}
    
    # Shapiro-Wilk Test (most powerful for small-medium samples)
    if n <= 5000:
        sw_stat, sw_p = stats.shapiro(clean_data)
        results['shapiro_wilk'] = {
            'statistic': sw_stat,
            'p_value': sw_p,
            'normal': sw_p > 0.05
        }
    
    # Anderson-Darling Test
    ad_result = stats.anderson(clean_data, dist='norm')
    results['anderson_darling'] = {
        'statistic': ad_result.statistic,
        'critical_values': ad_result.critical_values,
        'significance_levels': ad_result.significance_level,
        'normal': ad_result.statistic < ad_result.critical_values[2]  # 5% level
    }
    
    # Kolmogorov-Smirnov Test
    ks_stat, ks_p = stats.kstest(clean_data, 'norm', 
                                  args=(np.mean(clean_data), np.std(clean_data)))
    results['kolmogorov_smirnov'] = {
        'statistic': ks_stat,
        'p_value': ks_p,
        'normal': ks_p > 0.05
    }
    
    # D'Agostino-Pearson Test
    if n >= 8:
        dp_stat, dp_p = stats.normaltest(clean_data)
        results['dagostino_pearson'] = {
            'statistic': dp_stat,
            'p_value': dp_p,
            'normal': dp_p > 0.05
        }
    
    # Jarque-Bera Test
    jb_stat, jb_p = stats.jarque_bera(clean_data)
    results['jarque_bera'] = {
        'statistic': jb_stat,
        'p_value': jb_p,
        'normal': jb_p > 0.05
    }
    
    return results


def perform_comparison_tests(lgm_data, holocene_data, var_name):
    """Perform comprehensive non-parametric comparison tests"""
    results = {}
    
    # Clean data
    lgm_clean = lgm_data[~np.isnan(lgm_data)]
    holocene_clean = holocene_data[~np.isnan(holocene_data)]
    
    # Mann-Whitney U Test (primary non-parametric test)
    mw_stat, mw_p = stats.mannwhitneyu(lgm_clean, holocene_clean, 
                                        alternative='two-sided')
    results['mann_whitney_u'] = {
        'statistic': mw_stat,
        'p_value': mw_p,
        'significant': mw_p < 0.05
    }
    
    # Kruskal-Wallis H Test (alternative non-parametric)
    kw_stat, kw_p = stats.kruskal(lgm_clean, holocene_clean)
    results['kruskal_wallis'] = {
        'statistic': kw_stat,
        'p_value': kw_p,
        'significant': kw_p < 0.05
    }
    
    # Mood's Median Test
    mood_stat, mood_p, mood_med, mood_table = stats.median_test(lgm_clean, 
                                                                  holocene_clean)
    results['moods_median'] = {
        'statistic': mood_stat,
        'p_value': mood_p,
        'significant': mood_p < 0.05,
        'median': mood_med
    }
    
    # Wilcoxon Rank-Sum (equivalent to Mann-Whitney)
    wilcox_stat, wilcox_p = stats.ranksums(lgm_clean, holocene_clean)
    results['wilcoxon_ranksum'] = {
        'statistic': wilcox_stat,
        'p_value': wilcox_p,
        'significant': wilcox_p < 0.05
    }
    
    # Cliff's Delta (effect size)
    delta, delta_interp = cliffs_delta(lgm_clean, holocene_clean)
    results['cliffs_delta'] = {
        'delta': delta,
        'interpretation': delta_interp,
        'direction': 'LGM > Holocene' if delta > 0 else 'LGM < Holocene'
    }
    
    # Kolmogorov-Smirnov 2-sample test
    ks2_stat, ks2_p = stats.ks_2samp(lgm_clean, holocene_clean)
    results['ks_2sample'] = {
        'statistic': ks2_stat,
        'p_value': ks2_p,
        'significant': ks2_p < 0.05
    }
    
    # Epps-Singleton Test (for different distributions)
    es_stat, es_p = stats.epps_singleton_2samp(lgm_clean, holocene_clean)
    results['epps_singleton'] = {
        'statistic': es_stat,
        'p_value': es_p,
        'significant': es_p < 0.05
    }
    
    return results


def perform_complexity_analysis(lgm_data, holocene_data, var_name):
    """Analyze complexity differences between groups"""
    results = {}
    
    # Clean data
    lgm_clean = lgm_data[~np.isnan(lgm_data)]
    holocene_clean = holocene_data[~np.isnan(holocene_data)]
    
    # Shannon Entropy (histogram-based)
    lgm_hist, _ = np.histogram(lgm_clean, bins=50, density=True)
    holocene_hist, _ = np.histogram(holocene_clean, bins=50, density=True)
    
    # Normalize
    lgm_hist = lgm_hist / lgm_hist.sum()
    holocene_hist = holocene_hist / holocene_hist.sum()
    
    results['shannon_entropy'] = {
        'lgm': entropy(lgm_hist + 1e-10),
        'holocene': entropy(holocene_hist + 1e-10),
        'difference': entropy(lgm_hist + 1e-10) - entropy(holocene_hist + 1e-10)
    }
    
    # Approximate Entropy
    try:
        lgm_apen = approximate_entropy(lgm_clean)
        holocene_apen = approximate_entropy(holocene_clean)
        results['approximate_entropy'] = {
            'lgm': lgm_apen,
            'holocene': holocene_apen,
            'difference': lgm_apen - holocene_apen
        }
    except:
        results['approximate_entropy'] = {'error': 'Calculation failed'}
    
    # Sample Entropy
    try:
        lgm_sampen = sample_entropy(lgm_clean)
        holocene_sampen = sample_entropy(holocene_clean)
        results['sample_entropy'] = {
            'lgm': lgm_sampen,
            'holocene': holocene_sampen,
            'difference': lgm_sampen - holocene_sampen
        }
    except:
        results['sample_entropy'] = {'error': 'Calculation failed'}
    
    # Coefficient of Variation (normalized variability)
    results['coefficient_variation'] = {
        'lgm': np.std(lgm_clean) / np.mean(lgm_clean) if np.mean(lgm_clean) != 0 else np.nan,
        'holocene': np.std(holocene_clean) / np.mean(holocene_clean) if np.mean(holocene_clean) != 0 else np.nan
    }
    
    # Interquartile Range / Median (robust measure of spread)
    results['robust_spread'] = {
        'lgm': stats.iqr(lgm_clean) / np.median(lgm_clean) if np.median(lgm_clean) != 0 else np.nan,
        'holocene': stats.iqr(holocene_clean) / np.median(holocene_clean) if np.median(holocene_clean) != 0 else np.nan
    }
    
    # Range normalized by median
    results['normalized_range'] = {
        'lgm': (np.max(lgm_clean) - np.min(lgm_clean)) / np.median(lgm_clean) if np.median(lgm_clean) != 0 else np.nan,
        'holocene': (np.max(holocene_clean) - np.min(holocene_clean)) / np.median(holocene_clean) if np.median(holocene_clean) != 0 else np.nan
    }
    
    return results


def generate_interpretation(var_name, normality_lgm, normality_holocene, 
                           comparison, complexity):
    """Generate human-readable interpretation"""
    interpretation = []
    
    interpretation.append(f"{'='*80}")
    interpretation.append(f"VARIABLE: {var_name.upper()}")
    interpretation.append(f"{'='*80}\n")
    
    # Normality Assessment
    interpretation.append("NORMALITY ASSESSMENT:")
    interpretation.append("-" * 60)
    
    for group, norm_results in [("LGM", normality_lgm), 
                                 ("Late Holocene", normality_holocene)]:
        interpretation.append(f"\n{group}:")
        if 'error' in norm_results:
            interpretation.append(f"  ERROR: {norm_results['error']}")
            continue
            
        normal_count = sum(1 for test in norm_results.values() 
                          if isinstance(test, dict) and test.get('normal', False))
        total_tests = len([t for t in norm_results.values() if isinstance(t, dict)])
        
        interpretation.append(f"  Tests indicating normality: {normal_count}/{total_tests}")
        
        for test_name, test_result in norm_results.items():
            if isinstance(test_result, dict):
                if 'p_value' in test_result:
                    interpretation.append(
                        f"  - {test_name}: p={test_result['p_value']:.4f} "
                        f"({'NORMAL' if test_result['normal'] else 'NON-NORMAL'})"
                    )
                elif test_name == 'anderson_darling':
                    interpretation.append(
                        f"  - {test_name}: stat={test_result['statistic']:.4f} "
                        f"({'NORMAL' if test_result['normal'] else 'NON-NORMAL'})"
                    )
        
        interpretation.append(f"  CONCLUSION: Data is {'NORMAL' if normal_count > total_tests/2 else 'NON-NORMAL'}")
    
    # Comparison Tests
    interpretation.append("\n" + "="*60)
    interpretation.append("GROUP COMPARISON (NON-PARAMETRIC TESTS):")
    interpretation.append("-" * 60 + "\n")
    
    sig_count = sum(1 for test in comparison.values() 
                   if isinstance(test, dict) and test.get('significant', False))
    total_comp_tests = len([t for t in comparison.values() 
                           if isinstance(t, dict) and 'significant' in t])
    
    interpretation.append(f"Tests showing significant difference: {sig_count}/{total_comp_tests}\n")
    
    for test_name, test_result in comparison.items():
        if isinstance(test_result, dict) and 'p_value' in test_result:
            interpretation.append(
                f"{test_name}: p={test_result['p_value']:.6f} "
                f"({'SIGNIFICANT' if test_result['significant'] else 'NOT SIGNIFICANT'})"
            )
    
    # Effect Size
    if 'cliffs_delta' in comparison:
        cd = comparison['cliffs_delta']
        interpretation.append(f"\nEffect Size (Cliff's Delta):")
        interpretation.append(f"  δ = {cd['delta']:.4f} ({cd['interpretation']} effect)")
        interpretation.append(f"  Direction: {cd['direction']}")
    
    # Overall Conclusion
    interpretation.append("\n" + "="*60)
    interpretation.append("STATISTICAL CONCLUSION:")
    interpretation.append("-" * 60)
    if sig_count >= total_comp_tests * 0.7:  # 70% of tests
        interpretation.append(
            f"✓ STRONG EVIDENCE: LGM and Late Holocene {var_name} distributions "
            f"are SIGNIFICANTLY DIFFERENT ({sig_count}/{total_comp_tests} tests)"
        )
    elif sig_count >= total_comp_tests * 0.5:
        interpretation.append(
            f"~ MODERATE EVIDENCE: LGM and Late Holocene {var_name} distributions "
            f"show differences ({sig_count}/{total_comp_tests} tests)"
        )
    else:
        interpretation.append(
            f"✗ WEAK EVIDENCE: LGM and Late Holocene {var_name} distributions "
            f"are NOT significantly different ({sig_count}/{total_comp_tests} tests)"
        )
    
    # Complexity Analysis
    interpretation.append("\n" + "="*60)
    interpretation.append("COMPLEXITY ANALYSIS:")
    interpretation.append("-" * 60 + "\n")
    
    if 'shannon_entropy' in complexity:
        se = complexity['shannon_entropy']
        interpretation.append(f"Shannon Entropy:")
        interpretation.append(f"  LGM: {se['lgm']:.4f}")
        interpretation.append(f"  Holocene: {se['holocene']:.4f}")
        interpretation.append(f"  Difference: {se['difference']:.4f}")
        interpretation.append(f"  {'LGM is MORE complex' if se['difference'] > 0 else 'Holocene is MORE complex'}\n")
    
    if 'approximate_entropy' in complexity and 'error' not in complexity['approximate_entropy']:
        ae = complexity['approximate_entropy']
        interpretation.append(f"Approximate Entropy (regularity):")
        interpretation.append(f"  LGM: {ae['lgm']:.4f}")
        interpretation.append(f"  Holocene: {ae['holocene']:.4f}")
        interpretation.append(f"  {'LGM is MORE regular' if ae['difference'] < 0 else 'Holocene is MORE regular'}\n")
    
    if 'coefficient_variation' in complexity:
        cv = complexity['coefficient_variation']
        interpretation.append(f"Coefficient of Variation (relative variability):")
        interpretation.append(f"  LGM: {cv['lgm']:.4f}")
        interpretation.append(f"  Holocene: {cv['holocene']:.4f}\n")
    
    interpretation.append("\n" + "="*80 + "\n\n")
    
    return "\n".join(interpretation)


def main():
    # Set up paths
    data_dir = Path("../processed_data")
    stats_dir = Path("../stats")
    stats_dir.mkdir(parents=True, exist_ok=True)
    
    # Load data
    print("Loading data...")
    lgm_df = pd.read_csv(data_dir / "lgm_data.csv")
    holocene_df = pd.read_csv(data_dir / "late_holocene_data.csv")
    
    print(f"LGM data: {len(lgm_df)} records")
    print(f"Late Holocene data: {len(holocene_df)} records\n")
    
    # Variables to analyze
    variables = ['sst', 'sss', 'd18osw']
    
    # Store all results
    all_interpretations = []
    
    # Analyze each variable
    for var in variables:
        print(f"Analyzing {var}...")
        
        lgm_var = lgm_df[var].values
        holocene_var = holocene_df[var].values
        
        # Perform tests
        norm_lgm = perform_normality_tests(lgm_var, var, "LGM")
        norm_holocene = perform_normality_tests(holocene_var, var, "Late Holocene")
        comparison = perform_comparison_tests(lgm_var, holocene_var, var)
        complexity = perform_complexity_analysis(lgm_var, holocene_var, var)
        
        # Generate interpretation
        interpretation = generate_interpretation(
            var, norm_lgm, norm_holocene, comparison, complexity
        )
        
        all_interpretations.append(interpretation)
        
        # Save individual file
        output_file = stats_dir / f"{var}_statistical_analysis.txt"
        with open(output_file, 'w') as f:
            f.write(interpretation)
        
        print(f"  ✓ Saved to {output_file}")
    
    # Save combined report
    combined_file = stats_dir / "complete_statistical_test_report.txt"
    with open(combined_file, 'w') as f:
        f.write("RIGOROUS STATISTICAL ANALYSIS\n")
        f.write("LGM vs Late Holocene Comparison\n")
        f.write("="*80 + "\n\n")
        f.write("This analysis includes:\n")
        f.write("- Multiple normality tests (Shapiro-Wilk, Anderson-Darling, K-S, etc.)\n")
        f.write("- Non-parametric comparisons (Mann-Whitney U, Kruskal-Wallis, etc.)\n")
        f.write("- Effect sizes (Cliff's Delta)\n")
        f.write("- Complexity measures (entropy, regularity)\n\n")
        f.write("="*80 + "\n\n")
        f.write("\n\n".join(all_interpretations))
    
    print(f"\n✓ Complete report saved to {combined_file}")
    print("\nAnalysis complete!")


if __name__ == "__main__":
    main()
