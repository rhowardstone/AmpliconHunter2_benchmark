#!/usr/bin/env python3
"""
Generate publication-quality figures from AmpliconHunter benchmarking results.
Version 4 with final adjustments: base-2 log scales, fixed legends, and proper layouts.
"""

import os
import sys
import re
import glob
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D
import matplotlib.ticker as ticker
import seaborn as sns
from scipy import stats
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# Set style for publication - increase font sizes slightly
plt.style.use('seaborn-v0_8-whitegrid')
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.size'] = 12
plt.rcParams['axes.labelsize'] = 13
plt.rcParams['axes.labelweight'] = 'bold'
plt.rcParams['axes.titlesize'] = 14
plt.rcParams['axes.titleweight'] = 'bold'
plt.rcParams['xtick.labelsize'] = 11
plt.rcParams['ytick.labelsize'] = 11
plt.rcParams['legend.fontsize'] = 11
plt.rcParams['figure.titlesize'] = 16
plt.rcParams['figure.titleweight'] = 'bold'
plt.rcParams['lines.linewidth'] = 2.0
plt.rcParams['lines.markersize'] = 0  # No markers
plt.rcParams['axes.linewidth'] = 1.5
plt.rcParams['axes.grid'] = True
plt.rcParams['grid.alpha'] = 0.3
plt.rcParams['grid.linestyle'] = '--'

# Define consistent colors for all versions
# IMPORTANT: Internal naming vs Display naming
# Internal (data files): AHv1, AHv2, AHv3, AHv4
# Display (plots):        AHv1.1, AHv2.α, AHv2.β, AHv2
VERSION_COLORS = {
    'AHv1': '#2E86AB',   # Nice blue     -> displays as AHv1.1
    'AHv2': '#A23B72',   # Purple/magenta -> displays as AHv2.α
    'AHv3': '#F18F01',   # Orange        -> displays as AHv2.β
    'AHv4': '#C73E1D'    # Red           -> displays as AHv2
}

# Display name mapping for legends and labels
DISPLAY_NAMES = {
    'AHv1': 'AHv1.1',
    'AHv2': 'AHv2.α',
    'AHv3': 'AHv2.β', 
    'AHv4': 'AHv2'
}

def extract_metrics_from_log(log_file):
    """Extract comprehensive metrics from a log file"""
    metrics = {}
    
    if not os.path.exists(log_file):
        return metrics
    
    with open(log_file, 'r') as f:
        content = f.read()
    
    # Extract various metrics using regex
    patterns = {
        'user_time': r'User time \(seconds\): ([\d.]+)',
        'system_time': r'System time \(seconds\): ([\d.]+)',
        'elapsed_time': r'Elapsed \(wall clock\) time.*: ([\d:]+\.?\d*)',
        'cpu_percent': r'Percent of CPU this job got: (\d+)%',
        'max_rss_kb': r'Maximum resident set size \(kbytes\): (\d+)',
        'major_faults': r'Major \(requiring I/O\) page faults: (\d+)',
        'minor_faults': r'Minor \(reclaiming a frame\) page faults: (\d+)',
        'voluntary_switches': r'Voluntary context switches: (\d+)',
        'involuntary_switches': r'Involuntary context switches: (\d+)',
        'fs_inputs': r'File system inputs: (\d+)',
        'fs_outputs': r'File system outputs: (\d+)',
        'swaps': r'Swaps: (\d+)',
    }
    
    for key, pattern in patterns.items():
        match = re.search(pattern, content)
        if match:
            if key == 'elapsed_time':
                metrics[key] = parse_time_to_seconds(match.group(1))
            else:
                value = match.group(1)
                metrics[key] = float(value) if '.' in value else int(value)
    
    return metrics

def parse_time_to_seconds(time_str):
    """Convert time format (h:mm:ss or m:ss.ss) to seconds"""
    parts = time_str.split(':')
    if len(parts) == 3:  # h:mm:ss
        return int(parts[0]) * 3600 + int(parts[1]) * 60 + float(parts[2])
    elif len(parts) == 2:  # m:ss.ss
        return int(parts[0]) * 60 + float(parts[1])
    else:
        return float(time_str)

def read_metric_file(filepath):
    """Read a metric value from a file"""
    try:
        with open(filepath, 'r') as f:
            value = f.read().strip()
            if value and value != 'TIMEOUT_OR_ERROR':
                # Handle time format
                if ':' in value:
                    return parse_time_to_seconds(value)
                return float(value)
    except:
        pass
    return None

def collect_test_data(results_dir, test_num):
    """Collect all data for a specific test including rich metrics"""
    data = []
    
    if test_num == 1:
        # Test 1: Input Size Scaling
        datasets = ['G006400', 'G012800', 'G025600', 'G051200', 'G102400', 'G204800']
        genome_counts = [6400, 12800, 25600, 51200, 102400, 204800]
        
        for version in ['AHv1', 'AHv2', 'AHv3', 'AHv4']:
            for dataset, genome_count in zip(datasets, genome_counts):
                for repeat in range(1, 6):
                    base_path = f"{results_dir}/test1/{version}/{dataset}_repeat{repeat}"
                    
                    real_time = read_metric_file(f"{base_path}.real_time")
                    max_memory = read_metric_file(f"{base_path}.max_memory")
                    log_metrics = extract_metrics_from_log(f"{base_path}.log")
                    
                    if real_time is not None:
                        row = {
                            'version': version,
                            'dataset': dataset,
                            'genome_count': genome_count,
                            'repeat': repeat,
                            'real_time': real_time,
                            'max_memory_kb': max_memory if max_memory else log_metrics.get('max_rss_kb'),
                            **log_metrics
                        }
                        data.append(row)
    
    elif test_num == 2:
        # Test 2: N-base experiment
        n_values = [0, 2, 4, 6]
        
        for version in ['AHv1', 'AHv2', 'AHv3', 'AHv4']:
            for n in n_values:
                for repeat in range(1, 6):
                    base_path = f"{results_dir}/test2/{version}/{n}N_repeat{repeat}"
                    
                    real_time = read_metric_file(f"{base_path}.real_time")
                    max_memory = read_metric_file(f"{base_path}.max_memory")
                    log_metrics = extract_metrics_from_log(f"{base_path}.log")
                    
                    if real_time is not None:
                        row = {
                            'version': version,
                            'n_bases': n,
                            'repeat': repeat,
                            'real_time': real_time,
                            'max_memory_kb': max_memory if max_memory else log_metrics.get('max_rss_kb'),
                            **log_metrics
                        }
                        data.append(row)
    
    elif test_num == 3:
        # Test 3: Mismatch variation
        mismatches = [0, 1, 2, 3, 4, 5, 6]
        
        for version in ['AHv1', 'AHv2', 'AHv3', 'AHv4']:
            for mm in mismatches:
                for repeat in range(1, 6):
                    base_path = f"{results_dir}/test3/{version}/mm{mm}_repeat{repeat}"
                    
                    real_time = read_metric_file(f"{base_path}.real_time")
                    max_memory = read_metric_file(f"{base_path}.max_memory")
                    log_metrics = extract_metrics_from_log(f"{base_path}.log")
                    
                    if real_time is not None:
                        row = {
                            'version': version,
                            'mismatches': mm,
                            'repeat': repeat,
                            'real_time': real_time,
                            'max_memory_kb': max_memory if max_memory else log_metrics.get('max_rss_kb'),
                            **log_metrics
                        }
                        data.append(row)
    
    elif test_num == 4:
        # Test 4: Thread scaling
        thread_counts = [1, 2, 4, 8, 16, 32, 64, 96, 128, 160, 190]
        
        for version in ['AHv1', 'AHv2', 'AHv3', 'AHv4']:
            for threads in thread_counts:
                for repeat in range(1, 4):  # Only 3 repeats for test 4
                    base_path = f"{results_dir}/test4/{version}/t{threads}_repeat{repeat}"
                    
                    real_time = read_metric_file(f"{base_path}.real_time")
                    max_memory = read_metric_file(f"{base_path}.max_memory")
                    log_metrics = extract_metrics_from_log(f"{base_path}.log")
                    
                    if real_time is not None:
                        row = {
                            'version': version,
                            'threads': threads,
                            'repeat': repeat,
                            'real_time': real_time,
                            'max_memory_kb': max_memory if max_memory else log_metrics.get('max_rss_kb'),
                            **log_metrics
                        }
                        data.append(row)
    
    elif test_num == 5:
        # Test 5: Hot vs Cold Cache
        for version in ['AHv1', 'AHv2', 'AHv3', 'AHv4']:
            # Cold cache runs
            for repeat in range(1, 6):
                base_path = f"{results_dir}/test5/{version}/cold_repeat{repeat}"
                
                real_time = read_metric_file(f"{base_path}.real_time")
                max_memory = read_metric_file(f"{base_path}.max_memory")
                log_metrics = extract_metrics_from_log(f"{base_path}.log")
                
                if real_time is not None:
                    row = {
                        'version': version,
                        'cache_state': 'cold',
                        'repeat': repeat,
                        'real_time': real_time,
                        'max_memory_kb': max_memory if max_memory else log_metrics.get('max_rss_kb'),
                        **log_metrics
                    }
                    data.append(row)
            
            # Hot cache runs
            for repeat in range(1, 6):
                base_path = f"{results_dir}/test5/{version}/hot_repeat{repeat}"
                
                real_time = read_metric_file(f"{base_path}.real_time")
                max_memory = read_metric_file(f"{base_path}.max_memory")
                log_metrics = extract_metrics_from_log(f"{base_path}.log")
                
                if real_time is not None:
                    row = {
                        'version': version,
                        'cache_state': 'hot',
                        'repeat': repeat,
                        'real_time': real_time,
                        'max_memory_kb': max_memory if max_memory else log_metrics.get('max_rss_kb'),
                        **log_metrics
                    }
                    data.append(row)
    
    elif test_num == 6:
        # Test 6: Primer Pair Comparison
        primer_pairs = ['V3V4', 'Titan', 'V1V9']
        
        for version in ['AHv1', 'AHv2', 'AHv3', 'AHv4']:
            for primer in primer_pairs:
                for repeat in range(1, 6):
                    base_path = f"{results_dir}/test6/{version}/{primer}_repeat{repeat}"
                    
                    real_time = read_metric_file(f"{base_path}.real_time")
                    max_memory = read_metric_file(f"{base_path}.max_memory")
                    log_metrics = extract_metrics_from_log(f"{base_path}.log")
                    
                    if real_time is not None:
                        row = {
                            'version': version,
                            'primer_pair': primer,
                            'repeat': repeat,
                            'real_time': real_time,
                            'max_memory_kb': max_memory if max_memory else log_metrics.get('max_rss_kb'),
                            **log_metrics
                        }
                        data.append(row)
    
    return pd.DataFrame(data)

def calculate_ci95(data):
    """Calculate mean and 95% confidence interval"""
    n = len(data)
    if n == 0:
        return np.nan, np.nan, np.nan
    mean = np.mean(data)
    if n == 1:
        return mean, mean, mean
    sem = stats.sem(data)
    ci = sem * stats.t.ppf(0.975, n-1)
    return mean, mean-ci, mean+ci

def plot_with_ci(ax, x_values, y_values_list, label, color, linestyle='-', alpha_ci=0.25):
    """Plot line with 95% confidence interval shading - NO MARKERS"""
    means = []
    lower_bounds = []
    upper_bounds = []
    
    for values in y_values_list:
        mean, lower, upper = calculate_ci95(values)
        means.append(mean)
        lower_bounds.append(lower)
        upper_bounds.append(upper)
    
    # Plot line and shaded CI - no markers
    line = ax.plot(x_values, means, color=color, label=label, 
                   linestyle=linestyle, linewidth=2.5)
    ax.fill_between(x_values, lower_bounds, upper_bounds, color=color, alpha=alpha_ci)
    
    return means, lower_bounds, upper_bounds, line[0]

def create_main_figure(output_dir, test_data):
    """Create main Figure 1: Runtime comparison for all 6 tests with unified legend"""
    fig = plt.figure(figsize=(16, 14))
    
    # Create gridspec with moderate spacing
    gs = gridspec.GridSpec(3, 3, figure=fig, hspace=0.3, wspace=0.25, 
                          width_ratios=[1, 1, 0.15])
    
    # Store all line objects for legend
    legend_lines = []
    legend_labels = []
    
    # Get unique versions from the data
    all_versions = ['AHv1', 'AHv2', 'AHv3', 'AHv4']
    
    # Test 1: Input Size Scaling
    ax = fig.add_subplot(gs[0, 0])
    df1 = test_data[1]
    
    for version in all_versions:
        version_df = df1[df1['version'] == version]
        if not version_df.empty:
            x_values = sorted(version_df['genome_count'].unique())
            y_values = [version_df[version_df['genome_count'] == x]['real_time'].values for x in x_values]
            
            _, _, _, line = plot_with_ci(ax, x_values, y_values, DISPLAY_NAMES[version], 
                                         VERSION_COLORS[version])
            if len(legend_lines) < 4:  # Only collect once
                legend_lines.append(line)
                legend_labels.append(DISPLAY_NAMES[version])
    
    ax.set_xlabel('Number of Genomes', fontweight='bold')
    ax.set_ylabel('Runtime (s)', fontweight='bold')
    ax.set_title('(A) Input Size Scaling', fontweight='bold', pad=10)
    ax.set_xscale('log', base=2)
    ax.set_yscale('log', base=2)
    ax.grid(True, alpha=0.3, linestyle='--')
    # Set x-axis to show actual genome counts
    ax.set_xticks([6400, 12800, 25600, 51200, 102400, 204800])
    ax.set_xticklabels(['6400', '12800', '25600', '51200', '102400', '204800'])
    
    # Make tick labels bold
    for label in ax.get_xticklabels() + ax.get_yticklabels():
        label.set_fontweight('bold')
    
    # Test 2: N-base experiment
    ax = fig.add_subplot(gs[0, 1])
    df2 = test_data[2]
    
    for version in all_versions:
        version_df = df2[df2['version'] == version]
        if not version_df.empty:
            x_values = sorted(version_df['n_bases'].unique())
            y_values = [version_df[version_df['n_bases'] == x]['real_time'].values for x in x_values]
            
            plot_with_ci(ax, x_values, y_values, DISPLAY_NAMES[version], 
                        VERSION_COLORS[version])
    
    ax.set_xlabel('Number of N-bases', fontweight='bold')
    ax.set_ylabel('Runtime (s)', fontweight='bold')
    ax.set_title('(B) Primer Degeneracy Scaling', fontweight='bold', pad=10)
    ax.set_yscale('log', base=2)
    ax.grid(True, alpha=0.3, linestyle='--')
    ax.set_xticks([0, 2, 4, 6])
    
    for label in ax.get_xticklabels() + ax.get_yticklabels():
        label.set_fontweight('bold')
    
    # Test 3: Mismatch variation
    ax = fig.add_subplot(gs[1, 0])
    df3 = test_data[3]
    
    for version in all_versions:
        version_df = df3[df3['version'] == version]
        if not version_df.empty:
            x_values = sorted(version_df['mismatches'].unique())
            y_values = [version_df[version_df['mismatches'] == x]['real_time'].values for x in x_values]
            
            plot_with_ci(ax, x_values, y_values, DISPLAY_NAMES[version], 
                        VERSION_COLORS[version])
    
    ax.set_xlabel('Number of Mismatches', fontweight='bold')
    ax.set_ylabel('Runtime (s)', fontweight='bold')
    ax.set_title('(C) Mismatch Tolerance', fontweight='bold', pad=10)
    ax.set_yscale('log', base=2)
    ax.grid(True, alpha=0.3, linestyle='--')
    ax.set_xticks(range(7))
    
    for label in ax.get_xticklabels() + ax.get_yticklabels():
        label.set_fontweight('bold')
    
    # Test 4: Thread scaling
    ax = fig.add_subplot(gs[1, 1])
    df4 = test_data[4]
    
    for version in all_versions:
        version_df = df4[df4['version'] == version]
        if not version_df.empty:
            x_values = sorted(version_df['threads'].unique())
            y_values = [version_df[version_df['threads'] == x]['real_time'].values for x in x_values]
            
            plot_with_ci(ax, x_values, y_values, DISPLAY_NAMES[version], 
                        VERSION_COLORS[version])
    
    ax.set_xlabel('Number of Threads', fontweight='bold')
    ax.set_ylabel('Runtime (s)', fontweight='bold')
    ax.set_title('(D) Thread Scaling', fontweight='bold', pad=10)
    ax.set_xscale('log', base=2)
    ax.set_yscale('log', base=2)
    ax.grid(True, alpha=0.3, linestyle='--')
    ax.set_xticks([1, 2, 4, 8, 16, 32, 64, 128, 190])
    ax.set_xticklabels(['1', '2', '4', '8', '16', '32', '64', '128', '190'], fontweight='bold')
    
    for label in ax.get_yticklabels():
        label.set_fontweight('bold')
    
    # Test 5: Cache Performance - Group by cache state with error bars
    ax = fig.add_subplot(gs[2, 0])
    df5 = test_data[5]
    
    if not df5.empty:
        # Reorganize data to group by cache state
        cache_states = ['cold', 'hot']
        x = np.arange(len(cache_states))
        width = 0.18
        
        for i, version in enumerate(all_versions):
            version_means = []
            version_errors = []
            
            for state in cache_states:
                times = df5[(df5['version'] == version) & (df5['cache_state'] == state)]['real_time'].values
                if len(times) > 0:
                    mean, lower, upper = calculate_ci95(times)
                    version_means.append(mean)
                    version_errors.append(upper - mean)
                else:
                    version_means.append(0)
                    version_errors.append(0)
            
            if any(m > 0 for m in version_means):
                offset = (i - len(all_versions)/2 + 0.5) * width
                ax.bar(x + offset, version_means, width, yerr=version_errors,
                      color=VERSION_COLORS[version], capsize=4, alpha=0.8,
                      error_kw={'linewidth': 1.5, 'elinewidth': 1.5})
        
        ax.set_xlabel('Cache State', fontweight='bold')
        ax.set_ylabel('Runtime (s)', fontweight='bold')
        ax.set_title('(E) Cache Performance', fontweight='bold', pad=10)
        ax.set_xticks(x)
        ax.set_xticklabels(['Cold', 'Hot'], fontweight='bold')
        ax.grid(True, alpha=0.3, linestyle='--', axis='y')
        
        for label in ax.get_yticklabels():
            label.set_fontweight('bold')
    
    # Test 6: Primer Pair Comparison - Group by primer with error bars
    ax = fig.add_subplot(gs[2, 1])
    df6 = test_data[6]
    
    if not df6.empty:
        primers = sorted(df6['primer_pair'].unique())
        x = np.arange(len(primers))
        width = 0.18
        
        for i, version in enumerate(all_versions):
            version_means = []
            version_errors = []
            
            for primer in primers:
                times = df6[(df6['version'] == version) & (df6['primer_pair'] == primer)]['real_time'].values
                if len(times) > 0:
                    mean, lower, upper = calculate_ci95(times)
                    version_means.append(mean)
                    version_errors.append(upper - mean)
                else:
                    version_means.append(0)
                    version_errors.append(0)
            
            if any(m > 0 for m in version_means):
                offset = (i - len(all_versions)/2 + 0.5) * width
                ax.bar(x + offset, version_means, width, yerr=version_errors,
                      color=VERSION_COLORS[version], capsize=4, alpha=0.8,
                      error_kw={'linewidth': 1.5, 'elinewidth': 1.5})
        
        ax.set_xlabel('Primer Pair', fontweight='bold')
        ax.set_ylabel('Runtime (s)', fontweight='bold')
        ax.set_title('(F) Primer Pair Performance', fontweight='bold', pad=10)
        ax.set_xticks(x)
        ax.set_xticklabels(primers, fontweight='bold')
        ax.grid(True, alpha=0.3, linestyle='--', axis='y')
        
        for label in ax.get_yticklabels():
            label.set_fontweight('bold')
    
    # Add unified legend in the third column
    ax_legend = fig.add_subplot(gs[:, 2])
    ax_legend.axis('off')
    
    # Create legend with version names and colors
    legend = ax_legend.legend(legend_lines, legend_labels, 
                              loc='center', fontsize=12,
                              title='Version', title_fontsize=13,
                              frameon=True, fancybox=True, shadow=True)
    
    # Make legend text bold
    legend.get_title().set_fontweight('bold')
    plt.setp(legend.get_texts(), fontweight='bold')
    
    # Adjust layout
    plt.tight_layout()
    
    # Save in EPS format only
    fig.savefig(f'{output_dir}/Figure1_runtime_comparison.jpg', dpi=300, bbox_inches='tight')
    
    plt.close()
    print("Created Figure 1: Runtime comparison")

def create_memory_supplemental(output_dir, test_data):
    """Create Supplemental Figure S1: Memory usage analysis with 3x2 layout like Figure 1"""
    fig = plt.figure(figsize=(16, 14))
    
    # Create gridspec with same layout as Figure 1
    gs = gridspec.GridSpec(3, 3, figure=fig, hspace=0.3, wspace=0.25,
                          width_ratios=[1, 1, 0.15])
    
    all_versions = ['AHv1', 'AHv2', 'AHv3', 'AHv4']
    
    # Store legend elements
    legend_lines = []
    legend_labels = []
    
    # Test 1: Memory scaling with input size
    ax = fig.add_subplot(gs[0, 0])
    df1 = test_data[1]
    
    for version in all_versions:
        version_df = df1[df1['version'] == version]
        if not version_df.empty and 'max_memory_kb' in version_df.columns:
            version_df = version_df.copy()
            version_df['max_memory_gb'] = version_df['max_memory_kb'] / (1024 * 1024)
            x_values = sorted(version_df['genome_count'].unique())
            y_values = [version_df[version_df['genome_count'] == x]['max_memory_gb'].values for x in x_values]
            
            _, _, _, line = plot_with_ci(ax, x_values, y_values, DISPLAY_NAMES[version], 
                                         VERSION_COLORS[version])
            if len(legend_lines) < 4:
                legend_lines.append(line)
                legend_labels.append(DISPLAY_NAMES[version])
    
    ax.set_xlabel('Number of Genomes', fontweight='bold')
    ax.set_ylabel('Peak Memory (GB)', fontweight='bold')
    ax.set_title('(A) Input Size', fontweight='bold', pad=10)
    ax.set_xscale('log', base=2)
    ax.set_yscale('log', base=2)
    ax.grid(True, alpha=0.3)
    # Set x-axis to show actual genome counts
    ax.set_xticks([6400, 12800, 25600, 51200, 102400, 204800])
    ax.set_xticklabels(['6400', '12800', '25600', '51200', '102400', '204800'])
    
    for label in ax.get_xticklabels() + ax.get_yticklabels():
        label.set_fontweight('bold')
    
    # Test 2: Memory vs N-bases
    ax = fig.add_subplot(gs[0, 1])
    df2 = test_data[2]
    
    for version in all_versions:
        version_df = df2[df2['version'] == version]
        if not version_df.empty and 'max_memory_kb' in version_df.columns:
            version_df = version_df.copy()
            version_df['max_memory_gb'] = version_df['max_memory_kb'] / (1024 * 1024)
            x_values = sorted(version_df['n_bases'].unique())
            y_values = [version_df[version_df['n_bases'] == x]['max_memory_gb'].values for x in x_values]
            
            plot_with_ci(ax, x_values, y_values, DISPLAY_NAMES[version], 
                        VERSION_COLORS[version])
    
    ax.set_xlabel('Number of N-bases', fontweight='bold')
    ax.set_ylabel('Peak Memory (GB)', fontweight='bold')
    ax.set_title('(B) Primer Degeneracy', fontweight='bold', pad=10)
    ax.grid(True, alpha=0.3)
    ax.set_xticks([0, 2, 4, 6])
    
    for label in ax.get_xticklabels() + ax.get_yticklabels():
        label.set_fontweight('bold')
    
    # Test 3: Memory vs Mismatches
    ax = fig.add_subplot(gs[1, 0])
    df3 = test_data[3]
    
    for version in all_versions:
        version_df = df3[df3['version'] == version]
        if not version_df.empty and 'max_memory_kb' in version_df.columns:
            version_df = version_df.copy()
            version_df['max_memory_gb'] = version_df['max_memory_kb'] / (1024 * 1024)
            x_values = sorted(version_df['mismatches'].unique())
            y_values = [version_df[version_df['mismatches'] == x]['max_memory_gb'].values for x in x_values]
            
            plot_with_ci(ax, x_values, y_values, DISPLAY_NAMES[version], 
                        VERSION_COLORS[version])
    
    ax.set_xlabel('Number of Mismatches', fontweight='bold')
    ax.set_ylabel('Peak Memory (GB)', fontweight='bold')
    ax.set_title('(C) Mismatch Tolerance', fontweight='bold', pad=10)
    ax.grid(True, alpha=0.3)
    ax.set_xticks(range(7))
    
    for label in ax.get_xticklabels() + ax.get_yticklabels():
        label.set_fontweight('bold')
    
    # Test 4: Memory vs Thread count
    ax = fig.add_subplot(gs[1, 1])
    df4 = test_data[4]
    
    for version in all_versions:
        version_df = df4[df4['version'] == version]
        if not version_df.empty and 'max_memory_kb' in version_df.columns:
            version_df = version_df.copy()
            version_df['max_memory_gb'] = version_df['max_memory_kb'] / (1024 * 1024)
            x_values = sorted(version_df['threads'].unique())
            y_values = [version_df[version_df['threads'] == x]['max_memory_gb'].values for x in x_values]
            
            plot_with_ci(ax, x_values, y_values, DISPLAY_NAMES[version], 
                        VERSION_COLORS[version])
    
    ax.set_xlabel('Number of Threads', fontweight='bold')
    ax.set_ylabel('Peak Memory (GB)', fontweight='bold')
    ax.set_title('(D) Thread Scaling', fontweight='bold', pad=10)
    ax.set_xscale('log', base=2)
    ax.grid(True, alpha=0.3)
    ax.set_xticks([1, 2, 4, 8, 16, 32, 64, 128, 190])
    ax.set_xticklabels(['1', '2', '4', '8', '16', '32', '64', '128', '190'], fontweight='bold')
    
    for label in ax.get_yticklabels():
        label.set_fontweight('bold')
    
    # Test 5: Memory for cache states - grouped bar with error bars
    ax = fig.add_subplot(gs[2, 0])
    df5 = test_data[5]
    
    if not df5.empty and 'max_memory_kb' in df5.columns:
        cache_states = ['cold', 'hot']
        x = np.arange(len(cache_states))
        width = 0.18
        
        for i, version in enumerate(all_versions):
            version_means = []
            version_errors = []
            
            for state in cache_states:
                df5_filtered = df5[(df5['version'] == version) & (df5['cache_state'] == state)].copy()
                if not df5_filtered.empty:
                    df5_filtered['max_memory_gb'] = df5_filtered['max_memory_kb'] / (1024 * 1024)
                    mean, lower, upper = calculate_ci95(df5_filtered['max_memory_gb'].values)
                    version_means.append(mean)
                    version_errors.append(upper - mean)
                else:
                    version_means.append(0)
                    version_errors.append(0)
            
            if any(m > 0 for m in version_means):
                offset = (i - len(all_versions)/2 + 0.5) * width
                ax.bar(x + offset, version_means, width, yerr=version_errors,
                      color=VERSION_COLORS[version], alpha=0.8, capsize=4,
                      error_kw={'linewidth': 1.5, 'elinewidth': 1.5})
        
        ax.set_xlabel('Cache State', fontweight='bold')
        ax.set_ylabel('Peak Memory (GB)', fontweight='bold')
        ax.set_title('(E) Cache State', fontweight='bold', pad=10)
        ax.set_xticks(x)
        ax.set_xticklabels(['Cold', 'Hot'], fontweight='bold')
        ax.grid(True, alpha=0.3, axis='y')
        
        for label in ax.get_yticklabels():
            label.set_fontweight('bold')
    
    # Test 6: Memory for different primers - grouped bar with error bars
    ax = fig.add_subplot(gs[2, 1])
    df6 = test_data[6]
    
    if not df6.empty and 'max_memory_kb' in df6.columns:
        primers = sorted(df6['primer_pair'].unique())
        x = np.arange(len(primers))
        width = 0.18
        
        for i, version in enumerate(all_versions):
            version_means = []
            version_errors = []
            
            for primer in primers:
                df6_filtered = df6[(df6['version'] == version) & (df6['primer_pair'] == primer)].copy()
                if not df6_filtered.empty:
                    df6_filtered['max_memory_gb'] = df6_filtered['max_memory_kb'] / (1024 * 1024)
                    mean, lower, upper = calculate_ci95(df6_filtered['max_memory_gb'].values)
                    version_means.append(mean)
                    version_errors.append(upper - mean)
                else:
                    version_means.append(0)
                    version_errors.append(0)
            
            if any(m > 0 for m in version_means):
                offset = (i - len(all_versions)/2 + 0.5) * width
                ax.bar(x + offset, version_means, width, yerr=version_errors,
                      color=VERSION_COLORS[version], alpha=0.8, capsize=4,
                      error_kw={'linewidth': 1.5, 'elinewidth': 1.5})
        
        ax.set_xlabel('Primer Pair', fontweight='bold')
        ax.set_ylabel('Peak Memory (GB)', fontweight='bold')
        ax.set_title('(F) Primer Pairs', fontweight='bold', pad=10)
        ax.set_xticks(x)
        ax.set_xticklabels(primers, fontweight='bold')
        ax.grid(True, alpha=0.3, axis='y')
        
        for label in ax.get_yticklabels():
            label.set_fontweight('bold')
    
    # Add unified legend
    ax_legend = fig.add_subplot(gs[:, 2])
    ax_legend.axis('off')
    legend = ax_legend.legend(legend_lines, legend_labels, 
                              loc='center', fontsize=12,
                              title='Version', title_fontsize=13,
                              frameon=True, fancybox=True, shadow=True)
    legend.get_title().set_fontweight('bold')
    plt.setp(legend.get_texts(), fontweight='bold')
    
    plt.tight_layout()
    fig.savefig(f'{output_dir}/FigureS1_memory_analysis.jpg', dpi=300, bbox_inches='tight')
    plt.close()
    print("Created Figure S1: Memory analysis")

def create_efficiency_analysis(output_dir, test_data):
    """Create Figure S2: Efficiency analysis with proper y-axis limits and bold text"""
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    all_versions = ['AHv1', 'AHv2', 'AHv3', 'AHv4']
    
    # Input size scaling efficiency
    ax = axes[0, 0]
    df1 = test_data[1]
    
    max_efficiency = 0
    for version in all_versions:
        version_df = df1[df1['version'] == version]
        if not version_df.empty:
            x_values = sorted(version_df['genome_count'].unique())
            
            # Calculate average runtime for each genome count
            runtimes = []
            for x in x_values:
                mean_runtime = version_df[version_df['genome_count'] == x]['real_time'].mean()
                runtimes.append(mean_runtime)
            
            if len(runtimes) > 1 and runtimes[0] > 0:
                # Calculate scaling efficiency relative to smallest dataset
                baseline = runtimes[0]
                expected = [baseline * (x / x_values[0]) for x in x_values]
                efficiency = [expected[i] / runtimes[i] * 100 if runtimes[i] > 0 else 0 
                             for i in range(len(runtimes))]
                
                ax.plot(x_values, efficiency, color=VERSION_COLORS[version], label=DISPLAY_NAMES[version], 
                       linewidth=2.5)
                max_efficiency = max(max_efficiency, max(efficiency))
    
    ax.axhline(y=100, color='black', linestyle='--', alpha=0.5, label='Perfect Scaling')
    ax.set_xlabel('Number of Genomes', fontweight='bold')
    ax.set_ylabel('Scaling Efficiency (%)', fontweight='bold')
    ax.set_title('(A) Input Size Scaling Efficiency', fontweight='bold')
    ax.set_xscale('log', base=2)
    ax.grid(True, alpha=0.3)
    # Set x-axis to show actual genome counts
    ax.set_xticks([6400, 12800, 25600, 51200, 102400, 204800])
    ax.set_xticklabels(['6400', '12800', '25600', '51200', '102400', '204800'])
    legend = ax.legend(loc='best', fontsize=10, framealpha=0.95)
    plt.setp(legend.get_texts(), fontweight='bold')
    ax.set_ylim(0, min(150, max_efficiency * 1.1))  # Auto-fit with some headroom
    
    for label in ax.get_xticklabels() + ax.get_yticklabels():
        label.set_fontweight('bold')
    
    # Thread scaling efficiency
    ax = axes[0, 1]
    df4 = test_data[4]
    
    max_efficiency = 0
    for version in all_versions:
        version_df = df4[df4['version'] == version]
        if not version_df.empty:
            x_values = sorted(version_df['threads'].unique())
            
            # Calculate average runtime for each thread count
            runtimes = []
            for x in x_values:
                mean_runtime = version_df[version_df['threads'] == x]['real_time'].mean()
                runtimes.append(mean_runtime)
            
            # Calculate parallel efficiency
            if len(runtimes) > 0 and runtimes[0] > 0:
                speedup = [runtimes[0] / r if r > 0 else 0 for r in runtimes]
                efficiency = [speedup[i] / x_values[i] * 100 for i in range(len(speedup))]
                
                ax.plot(x_values, efficiency, color=VERSION_COLORS[version], label=DISPLAY_NAMES[version],
                       linewidth=2.5)
                max_efficiency = max(max_efficiency, max(efficiency))
    
    ax.axhline(y=100, color='black', linestyle='--', alpha=0.5, label='Perfect Scaling')
    ax.set_xlabel('Number of Threads', fontweight='bold')
    ax.set_ylabel('Parallel Efficiency (%)', fontweight='bold')
    ax.set_title('(B) Thread Scaling Efficiency', fontweight='bold')
    ax.set_xscale('log', base=2)
    ax.grid(True, alpha=0.3)
    legend = ax.legend(loc='best', fontsize=10, framealpha=0.95)
    plt.setp(legend.get_texts(), fontweight='bold')
    ax.set_ylim(0, min(150, max_efficiency * 1.1))
    
    for label in ax.get_xticklabels() + ax.get_yticklabels():
        label.set_fontweight('bold')
    
    # Speedup comparison
    ax = axes[1, 0]
    
    # Calculate speedup of each version relative to AHv1
    baseline_version = 'AHv1'
    if baseline_version in df1['version'].values:
        baseline_df = df1[df1['version'] == baseline_version]
        
        for version in all_versions:
            if version != baseline_version:
                version_df = df1[df1['version'] == version]
                if not version_df.empty:
                    x_values = sorted(set(version_df['genome_count'].unique()) & 
                                    set(baseline_df['genome_count'].unique()))
                    
                    speedups = []
                    for x in x_values:
                        baseline_time = baseline_df[baseline_df['genome_count'] == x]['real_time'].mean()
                        version_time = version_df[version_df['genome_count'] == x]['real_time'].mean()
                        if version_time > 0:
                            speedups.append(baseline_time / version_time)
                        else:
                            speedups.append(0)
                    
                    if speedups:
                        ax.plot(x_values, speedups, color=VERSION_COLORS[version], 
                               label=f'{DISPLAY_NAMES[version]} vs {DISPLAY_NAMES[baseline_version]}',
                               linewidth=2.5)
    
    ax.axhline(y=1, color='black', linestyle='--', alpha=0.5, label='No Speedup')
    ax.set_xlabel('Number of Genomes', fontweight='bold')
    ax.set_ylabel('Speedup Factor', fontweight='bold')
    ax.set_title('(C) Speedup Relative to AHv1', fontweight='bold')
    ax.set_xscale('log', base=2)
    ax.grid(True, alpha=0.3)
    # Set x-axis to show actual genome counts
    ax.set_xticks([6400, 12800, 25600, 51200, 102400, 204800])
    ax.set_xticklabels(['6400', '12800', '25600', '51200', '102400', '204800'])
    # Move legend inside plot area, on right side, slightly above center
    legend = ax.legend(loc='center right', bbox_to_anchor=(0.98, 0.55), 
                      bbox_transform=ax.transAxes, fontsize=10, framealpha=0.95)
    plt.setp(legend.get_texts(), fontweight='bold')
    
    for label in ax.get_xticklabels() + ax.get_yticklabels():
        label.set_fontweight('bold')
    
    # Cache efficiency
    ax = axes[1, 1]
    df5 = test_data[5]
    
    if not df5.empty:
        cache_speedup = []
        
        for version in all_versions:
            cold_times = df5[(df5['version'] == version) & (df5['cache_state'] == 'cold')]['real_time'].values
            hot_times = df5[(df5['version'] == version) & (df5['cache_state'] == 'hot')]['real_time'].values
            
            if len(cold_times) > 0 and len(hot_times) > 0:
                cold_mean = np.mean(cold_times)
                hot_mean = np.mean(hot_times)
                speedup = cold_mean / hot_mean if hot_mean > 0 else 1
                
                cache_speedup.append({
                    'version': version,
                    'speedup': speedup,
                    'cold': cold_mean,
                    'hot': hot_mean
                })
        
        if cache_speedup:
            x = np.arange(len(cache_speedup))
            speedups = [d['speedup'] for d in cache_speedup]
            speedup_errors = []
            
            # Calculate error bars for speedup ratios
            for version in all_versions:
                cold_times = df5[(df5['version'] == version) & (df5['cache_state'] == 'cold')]['real_time'].values
                hot_times = df5[(df5['version'] == version) & (df5['cache_state'] == 'hot')]['real_time'].values
                
                if len(cold_times) > 0 and len(hot_times) > 0:
                    # Bootstrap confidence interval for ratio
                    ratios = []
                    for _ in range(1000):
                        cold_sample = np.random.choice(cold_times, len(cold_times), replace=True)
                        hot_sample = np.random.choice(hot_times, len(hot_times), replace=True)
                        ratios.append(np.mean(cold_sample) / np.mean(hot_sample))
                    ci_lower, ci_upper = np.percentile(ratios, [2.5, 97.5])
                    speedup_errors.append((ci_upper - ci_lower) / 2)
                else:
                    speedup_errors.append(0)
            
            colors = [VERSION_COLORS[d['version']] for d in cache_speedup]
            
            bars = ax.bar(x, speedups, yerr=speedup_errors[:len(speedups)], 
                          color=colors, alpha=0.8, capsize=5,
                          error_kw={'linewidth': 1.5, 'elinewidth': 1.5})
            
            # NO text labels on bars - removed completely
            
            ax.axhline(y=1, color='black', linestyle='--', alpha=0.5, label='No Cache Benefit')
            ax.set_xlabel('Version', fontweight='bold')
            ax.set_ylabel('Cache Speedup Factor', fontweight='bold')
            ax.set_title('(D) Cache Performance Benefit', fontweight='bold')
            ax.set_xticks(x)
            ax.set_xticklabels([DISPLAY_NAMES[d['version']] for d in cache_speedup], fontweight='bold')
            ax.grid(True, alpha=0.3, axis='y')
            ax.set_ylim(0, max(speedups) * 1.1)  # Adjust y limit without label headroom
            
            for label in ax.get_yticklabels():
                label.set_fontweight('bold')
    
    plt.tight_layout()
    fig.savefig(f'{output_dir}/FigureS2_efficiency_analysis.jpg', dpi=300, bbox_inches='tight')
    plt.close()
    print("Created Figure S2: Efficiency analysis")

def create_io_system_metrics(output_dir, test_data):
    """Create Figure S3: I/O and System Metrics with 3x2 layout like other figures"""
    fig = plt.figure(figsize=(16, 14))
    
    # Create gridspec with same layout as Figure 1 (3x2 with legend column)
    gs = GridSpec(3, 3, figure=fig, hspace=0.3, wspace=0.25,
                  width_ratios=[1, 1, 0.15])
    
    all_versions = ['AHv1', 'AHv2', 'AHv3', 'AHv4']
    
    # Store legend elements
    legend_lines = []
    legend_labels = []
    
    # File system inputs - only plot if data exists and differs
    ax = fig.add_subplot(gs[0, 0])
    df1 = test_data[1]
    if 'fs_inputs' in df1.columns:
        has_data = False
        for version in all_versions:
            version_df = df1[df1['version'] == version]
            if not version_df.empty and version_df['fs_inputs'].notna().any():
                x_values = sorted(version_df['genome_count'].unique())
                # Convert to GB
                y_values = []
                for x in x_values:
                    vals = version_df[version_df['genome_count'] == x]['fs_inputs'].values / (1024**3)
                    vals = vals[~np.isnan(vals)]
                    if len(vals) > 0:
                        y_values.append(vals)
                    else:
                        y_values.append([])
                
                if any(len(y) > 0 for y in y_values):
                    _, _, _, line = plot_with_ci(ax, x_values, y_values, version, 
                                                 VERSION_COLORS[version])
                    has_data = True
                    if len(legend_lines) < 4:
                        legend_lines.append(line)
                        legend_labels.append(DISPLAY_NAMES[version])
        
        if has_data:
            ax.set_xlabel('Number of Genomes', fontweight='bold')
            ax.set_ylabel('Read Volume (GB)', fontweight='bold')
            ax.set_title('(A) I/O Read Volume', fontweight='bold')
            ax.set_xscale('log', base=2)
            ax.set_yscale('log', base=2)
            ax.grid(True, alpha=0.3)
            # Set x-axis to show actual genome counts
            ax.set_xticks([6400, 12800, 25600, 51200, 102400, 204800])
            ax.set_xticklabels(['6400', '12800', '25600', '51200', '102400', '204800'])
        else:
            ax.text(0.5, 0.5, 'No I/O data available', ha='center', va='center', 
                   fontsize=12, fontweight='bold', transform=ax.transAxes)
            ax.set_title('(A) I/O Read Volume', fontweight='bold')
        
        for label in ax.get_xticklabels() + ax.get_yticklabels():
            label.set_fontweight('bold')
    
    # Context switches
    ax = fig.add_subplot(gs[0, 1])
    if 'voluntary_switches' in df1.columns and 'involuntary_switches' in df1.columns:
        has_data = False
        for version in all_versions:
            version_df = df1[df1['version'] == version].copy()
            if not version_df.empty:
                version_df['total_switches'] = version_df['voluntary_switches'] + version_df['involuntary_switches']
                if version_df['total_switches'].notna().any():
                    x_values = sorted(version_df['genome_count'].unique())
                    y_values = []
                    for x in x_values:
                        vals = version_df[version_df['genome_count'] == x]['total_switches'].values / 1000
                        vals = vals[~np.isnan(vals)]
                        if len(vals) > 0:
                            y_values.append(vals)
                        else:
                            y_values.append([])
                    
                    if any(len(y) > 0 for y in y_values):
                        plot_with_ci(ax, x_values, y_values, version, 
                                    VERSION_COLORS[version])
                        has_data = True
        
        if has_data:
            ax.set_xlabel('Number of Genomes', fontweight='bold')
            ax.set_ylabel('Context Switches (×10³)', fontweight='bold')
            ax.set_title('(B) Context Switching', fontweight='bold')
            ax.set_xscale('log', base=2)
            ax.set_yscale('log', base=2)
            ax.grid(True, alpha=0.3)
            # Set x-axis to show actual genome counts
            ax.set_xticks([6400, 12800, 25600, 51200, 102400, 204800])
            ax.set_xticklabels(['6400', '12800', '25600', '51200', '102400', '204800'])
        else:
            ax.text(0.5, 0.5, 'No context switch data', ha='center', va='center',
                   fontsize=12, fontweight='bold', transform=ax.transAxes)
            ax.set_title('(B) Context Switching', fontweight='bold')
        
        for label in ax.get_xticklabels() + ax.get_yticklabels():
            label.set_fontweight('bold')
    
    # Page faults
    ax = fig.add_subplot(gs[1, 0])
    if 'major_faults' in df1.columns:
        has_data = False
        for version in all_versions:
            version_df = df1[df1['version'] == version]
            if not version_df.empty and version_df['major_faults'].notna().any():
                x_values = sorted(version_df['genome_count'].unique())
                y_values = []
                for x in x_values:
                    vals = version_df[version_df['genome_count'] == x]['major_faults'].values
                    vals = vals[~np.isnan(vals)]
                    if len(vals) > 0 and any(v > 0 for v in vals):
                        y_values.append(vals)
                    else:
                        y_values.append([])
                
                if any(len(y) > 0 for y in y_values):
                    plot_with_ci(ax, x_values, y_values, DISPLAY_NAMES[version], 
                                VERSION_COLORS[version])
                    has_data = True
        
        if has_data:
            ax.set_xlabel('Number of Genomes', fontweight='bold')
            ax.set_ylabel('Major Page Faults', fontweight='bold')
            ax.set_title('(C) I/O Page Faults', fontweight='bold')
            ax.set_xscale('log', base=2)
            ax.set_yscale('symlog')
            ax.grid(True, alpha=0.3)
            # Set x-axis to show actual genome counts
            ax.set_xticks([6400, 12800, 25600, 51200, 102400, 204800])
            ax.set_xticklabels(['6400', '12800', '25600', '51200', '102400', '204800'])
        else:
            ax.text(0.5, 0.5, 'No page fault data', ha='center', va='center',
                   fontsize=12, fontweight='bold', transform=ax.transAxes)
            ax.set_title('(C) I/O Page Faults', fontweight='bold')
        
        for label in ax.get_xticklabels() + ax.get_yticklabels():
            label.set_fontweight('bold')
    
    # CPU utilization by thread count
    ax = fig.add_subplot(gs[1, 1])
    df4 = test_data[4]
    if 'cpu_percent' in df4.columns:
        has_data = False
        for version in all_versions:
            version_df = df4[df4['version'] == version]
            if not version_df.empty and version_df['cpu_percent'].notna().any():
                x_values = sorted(version_df['threads'].unique())
                y_values = []
                for x in x_values:
                    vals = version_df[version_df['threads'] == x]['cpu_percent'].values
                    vals = vals[~np.isnan(vals)]
                    if len(vals) > 0:
                        y_values.append(vals)
                    else:
                        y_values.append([])
                
                if any(len(y) > 0 for y in y_values):
                    plot_with_ci(ax, x_values, y_values, DISPLAY_NAMES[version], 
                                VERSION_COLORS[version])
                    has_data = True
        
        if has_data:
            ax.set_xlabel('Number of Threads', fontweight='bold')
            ax.set_ylabel('CPU Utilization (%)', fontweight='bold')
            ax.set_title('(D) CPU Utilization vs Threads', fontweight='bold')
            ax.set_xscale('log', base=2)
            ax.grid(True, alpha=0.3)
        else:
            ax.text(0.5, 0.5, 'No CPU utilization data', ha='center', va='center',
                   fontsize=12, fontweight='bold', transform=ax.transAxes)
            ax.set_title('(D) CPU Utilization vs Threads', fontweight='bold')
        
        for label in ax.get_xticklabels() + ax.get_yticklabels():
            label.set_fontweight('bold')
    
    # System vs User time ratio
    ax = fig.add_subplot(gs[2, 0])
    if 'user_time' in df1.columns and 'system_time' in df1.columns:
        has_data = False
        for version in all_versions:
            version_df = df1[df1['version'] == version].copy()
            if not version_df.empty and version_df['user_time'].notna().any() and version_df['system_time'].notna().any():
                version_df['sys_user_ratio'] = version_df['system_time'] / (version_df['user_time'] + version_df['system_time']) * 100
                x_values = sorted(version_df['genome_count'].unique())
                y_values = []
                for x in x_values:
                    vals = version_df[version_df['genome_count'] == x]['sys_user_ratio'].values
                    vals = vals[~np.isnan(vals)]
                    if len(vals) > 0:
                        y_values.append(vals)
                    else:
                        y_values.append([])
                
                if any(len(y) > 0 for y in y_values):
                    plot_with_ci(ax, x_values, y_values, DISPLAY_NAMES[version], 
                                VERSION_COLORS[version])
                    has_data = True
        
        if has_data:
            ax.set_xlabel('Number of Genomes', fontweight='bold')
            ax.set_ylabel('System Time (%)', fontweight='bold')
            ax.set_title('(E) System vs User Time', fontweight='bold')
            ax.set_xscale('log', base=2)
            ax.grid(True, alpha=0.3)
            # Set x-axis to show actual genome counts
            ax.set_xticks([6400, 12800, 25600, 51200, 102400, 204800])
            ax.set_xticklabels(['6400', '12800', '25600', '51200', '102400', '204800'])
        else:
            ax.text(0.5, 0.5, 'No timing breakdown data', ha='center', va='center',
                   fontsize=12, fontweight='bold', transform=ax.transAxes)
            ax.set_title('(E) System vs User Time', fontweight='bold')
        
        for label in ax.get_xticklabels() + ax.get_yticklabels():
            label.set_fontweight('bold')
    
    # Minor page faults (memory reclaims)
    ax = fig.add_subplot(gs[2, 1])
    if 'minor_faults' in df1.columns:
        has_data = False
        for version in all_versions:
            version_df = df1[df1['version'] == version]
            if not version_df.empty and version_df['minor_faults'].notna().any():
                x_values = sorted(version_df['genome_count'].unique())
                y_values = []
                for x in x_values:
                    vals = version_df[version_df['genome_count'] == x]['minor_faults'].values / 1000000
                    vals = vals[~np.isnan(vals)]
                    if len(vals) > 0:
                        y_values.append(vals)
                    else:
                        y_values.append([])
                
                if any(len(y) > 0 for y in y_values):
                    plot_with_ci(ax, x_values, y_values, DISPLAY_NAMES[version], 
                                VERSION_COLORS[version])
                    has_data = True
        
        if has_data:
            ax.set_xlabel('Number of Genomes', fontweight='bold')
            ax.set_ylabel('Minor Page Faults (×10⁶)', fontweight='bold')
            ax.set_title('(F) Memory Reclaims', fontweight='bold')
            ax.set_xscale('log', base=2)
            ax.set_yscale('log', base=2)
            ax.grid(True, alpha=0.3)
            # Set x-axis to show actual genome counts
            ax.set_xticks([6400, 12800, 25600, 51200, 102400, 204800])
            ax.set_xticklabels(['6400', '12800', '25600', '51200', '102400', '204800'])
        else:
            ax.text(0.5, 0.5, 'No minor fault data', ha='center', va='center',
                   fontsize=12, fontweight='bold', transform=ax.transAxes)
            ax.set_title('(F) Memory Reclaims', fontweight='bold')
        
        for label in ax.get_xticklabels() + ax.get_yticklabels():
            label.set_fontweight('bold')
    
    # Add unified legend in the third column
    ax_legend = fig.add_subplot(gs[:, 2])
    ax_legend.axis('off')
    
    if legend_lines:
        legend = ax_legend.legend(legend_lines, legend_labels,
                                  loc='center', fontsize=12,
                                  title='Version', title_fontsize=13,
                                  frameon=True, fancybox=True, shadow=True)
        legend.get_title().set_fontweight('bold')
        plt.setp(legend.get_texts(), fontweight='bold')
    
    plt.tight_layout()
    fig.savefig(f'{output_dir}/FigureS3_io_system_metrics.jpg', dpi=300, bbox_inches='tight')
    plt.close()
    print("Created Figure S3: I/O and system metrics")

def main():
    if len(sys.argv) < 2:
        print("Usage: python generate_publication_figures_v4.py <results_dir> [output_dir]")
        sys.exit(1)
    
    results_dir = sys.argv[1]
    output_dir = sys.argv[2] if len(sys.argv) > 2 else os.path.join(results_dir, 'publication_plots_v4')
    
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    print(f"Loading data from: {results_dir}")
    print(f"Saving figures to: {output_dir}")
    print("-" * 50)
    
    # Collect data for all 6 tests
    test_data = {}
    for test_num in range(1, 7):
        print(f"Loading Test {test_num} data...")
        test_data[test_num] = collect_test_data(results_dir, test_num)
        print(f"  Found {len(test_data[test_num])} data points")
    
    print("-" * 50)
    
    # Generate figures
    create_main_figure(output_dir, test_data)
    create_memory_supplemental(output_dir, test_data)
    create_efficiency_analysis(output_dir, test_data)
    create_io_system_metrics(output_dir, test_data)
    
    print("-" * 50)
    print("All figures generated successfully!")
    print(f"\nOutput directory: {output_dir}")
    print("\nGenerated files (EPS format):")
    print("Main Figure:")
    print("  - Figure1_runtime_comparison.jpg")
    print("\nSupplemental Figures:")
    print("  - FigureS1_memory_analysis.jpg")
    print("  - FigureS2_efficiency_analysis.jpg")
    print("  - FigureS3_io_system_metrics.jpg")

if __name__ == "__main__":
    main()