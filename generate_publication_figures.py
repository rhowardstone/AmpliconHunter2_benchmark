#!/usr/bin/env python3
"""
Generate publication figures from AmpliconHunter benchmarking results.
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
from matplotlib.ticker import ScalarFormatter, FuncFormatter
import seaborn as sns
from scipy import stats
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# Set style for publication
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
VERSION_COLORS = {
    'AHv1': '#2E86AB',   # Nice blue     -> displays as AHv1.1 (updated python ref)
    'AHv2': '#A23B72',   # Purple/magenta -> displays as AHv2.α (initial C version - uses 2bit)
    'AHv3': '#F18F01',   # Orange        -> displays as AHv2.β (AVX2 parallelization)
    'AHv4': '#C73E1D',   # Red           -> displays as AHv2.γ (chunked 16MB reading)
    'AHv5': '#1B998B'    # Teal/green    -> displays as AHv2 (with melting temp)
}

# Display name mapping for legends and labels
DISPLAY_NAMES = {
    'AHv1': 'AHv1.1',
    'AHv2': 'AHv2.α',
    'AHv3': 'AHv2.β',
    'AHv4': 'AHv2.γ',
    'AHv5': 'AHv2'
}

def setup_log_axis_labels(ax, axis='y'):
    """Set up axis to show actual numbers instead of powers for log scales"""
    formatter = ScalarFormatter()
    formatter.set_scientific(False)
    formatter.set_useOffset(False)

    if axis == 'y':
        ax.yaxis.set_major_formatter(formatter)
        ax.yaxis.set_minor_formatter(ticker.NullFormatter())
    else:
        ax.xaxis.set_major_formatter(formatter)
        ax.xaxis.set_minor_formatter(ticker.NullFormatter())

def extract_metrics_from_log(log_file):
    """Extract comprehensive metrics from a log file"""
    metrics = {}

    if not os.path.exists(log_file):
        return metrics

    with open(log_file, 'r') as f:
        content = f.read()

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

        for version in ['AHv1', 'AHv2', 'AHv3', 'AHv4', 'AHv5']:
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

        for version in ['AHv1', 'AHv2', 'AHv3', 'AHv4', 'AHv5']:
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

        for version in ['AHv1', 'AHv2', 'AHv3', 'AHv4', 'AHv5']:
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

        for version in ['AHv1', 'AHv2', 'AHv3', 'AHv4', 'AHv5']:
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
        for version in ['AHv1', 'AHv2', 'AHv3', 'AHv4', 'AHv5']:
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

        for version in ['AHv1', 'AHv2', 'AHv3', 'AHv4', 'AHv5']:
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
    all_versions = ['AHv1', 'AHv2', 'AHv3', 'AHv4', 'AHv5']

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
            if len(legend_lines) < 5:  # Only collect once for all 5 versions
                legend_lines.append(line)
                legend_labels.append(DISPLAY_NAMES[version])

    ax.set_xlabel('Number of Genomes', fontweight='bold')
    ax.set_ylabel('Runtime (s)', fontweight='bold')
    ax.set_title('(A) Input Size Scaling', fontweight='bold', pad=10)
    ax.set_xscale('log', base=2)
    ax.set_yscale('log', base=2)
    setup_log_axis_labels(ax, 'y')
    ax.grid(True, alpha=0.3, linestyle='--')
    # Set x-axis to show genome counts
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
    setup_log_axis_labels(ax, 'y')
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
    setup_log_axis_labels(ax, 'y')
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
    setup_log_axis_labels(ax, 'y')
    ax.grid(True, alpha=0.3, linestyle='--')
    ax.set_xticks([1, 2, 4, 8, 16, 32, 64, 128, 190])
    ax.set_xticklabels(['1', '2', '4', '8', '16', '32', '64', '128', '190'], fontweight='bold')

    for label in ax.get_yticklabels():
        label.set_fontweight('bold')

    # Test 5: Buffering Effect - Group by buffering state with error bars
    ax = fig.add_subplot(gs[2, 0])
    df5 = test_data[5]

    if not df5.empty:
        # Reorganize data to group by buffering state
        buffer_states = ['cold', 'hot']
        x = np.arange(len(buffer_states))
        width = 0.15

        for i, version in enumerate(all_versions):
            version_means = []
            version_errors = []

            for state in buffer_states:
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

        ax.set_xlabel('File Buffering', fontweight='bold')
        ax.set_ylabel('Runtime (s)', fontweight='bold')
        ax.set_title('(E) File Buffering Effect', fontweight='bold', pad=10)
        ax.set_xticks(x)
        ax.set_xticklabels(['Unbuffered', 'Buffered'], fontweight='bold')
        ax.grid(True, alpha=0.3, linestyle='--', axis='y')

        for label in ax.get_yticklabels():
            label.set_fontweight('bold')

    # Test 6: Primer Pair Comparison - Grouped by primer with error bars
    ax = fig.add_subplot(gs[2, 1])
    df6 = test_data[6]

    if not df6.empty:
        primers = sorted(df6['primer_pair'].unique())
        x = np.arange(len(primers))
        width = 0.15

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

    # Save in JPG format
    fig.savefig(f'{output_dir}/Figure1_runtime_comparison.jpg', dpi=300, bbox_inches='tight')
    plt.close()
    print("Created Figure 1: Runtime comparison")

def create_figure1(output_dir, test_data):
    """Create Figure 1: 2x2 layout with S1-A, S2-A, S3-B, S3-C"""
    fig = plt.figure(figsize=(14, 12))

    # Create gridspec with 2x3 layout (2x2 plots + legend column)
    gs = gridspec.GridSpec(2, 3, figure=fig, hspace=0.3, wspace=0.25,
                          width_ratios=[1, 1, 0.15])

    all_versions = ['AHv1', 'AHv2', 'AHv3', 'AHv4', 'AHv5']

    # Store legend elements
    legend_lines = []
    legend_labels = []

    # Panel A:
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
            if len(legend_lines) < 5:
                legend_lines.append(line)
                legend_labels.append(DISPLAY_NAMES[version])

    ax.set_xlabel('Number of Genomes', fontweight='bold')
    ax.set_ylabel('Peak Memory (GB)', fontweight='bold')
    ax.set_title('(A) Memory: Input Size', fontweight='bold', pad=10)
    ax.set_xscale('log', base=2)
    ax.set_yscale('log', base=2)
    setup_log_axis_labels(ax, 'y')
    ax.grid(True, alpha=0.3)
    ax.set_xticks([6400, 12800, 25600, 51200, 102400, 204800])
    ax.set_xticklabels(['6400', '12800', '25600', '51200', '102400', '204800'])

    for label in ax.get_xticklabels() + ax.get_yticklabels():
        label.set_fontweight('bold')

    # Panel B:
    ax = fig.add_subplot(gs[0, 1])
    df1 = test_data[1]

    max_efficiency = 0
    for version in all_versions:
        version_df = df1[df1['version'] == version]
        if not version_df.empty:
            x_values = sorted(version_df['genome_count'].unique())

            runtimes = []
            for x in x_values:
                mean_runtime = version_df[version_df['genome_count'] == x]['real_time'].mean()
                runtimes.append(mean_runtime)

            if len(runtimes) > 1 and runtimes[0] > 0:
                baseline = runtimes[0]
                expected = [baseline * (x / x_values[0]) for x in x_values]
                efficiency = [expected[i] / runtimes[i] * 100 if runtimes[i] > 0 else 0
                             for i in range(len(runtimes))]

                ax.plot(x_values, efficiency, color=VERSION_COLORS[version],
                       label=DISPLAY_NAMES[version], linewidth=2.5)
                max_efficiency = max(max_efficiency, max(efficiency))

    ax.axhline(y=100, color='black', linestyle='--', alpha=0.5, label='Perfect Scaling')
    ax.set_xlabel('Number of Genomes', fontweight='bold')
    ax.set_ylabel('Scaling Efficiency (%)', fontweight='bold')
    ax.set_title('(B) Efficiency: Input Size Scaling', fontweight='bold', pad=10)
    ax.set_xscale('log', base=2)
    ax.grid(True, alpha=0.3)
    ax.set_xticks([6400, 12800, 25600, 51200, 102400, 204800])
    ax.set_xticklabels(['6400', '12800', '25600', '51200', '102400', '204800'])
    ax.set_ylim(0, min(150, max_efficiency * 1.1))

    for label in ax.get_xticklabels() + ax.get_yticklabels():
        label.set_fontweight('bold')

    # Panel C:
    ax = fig.add_subplot(gs[1, 0])
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
                        plot_with_ci(ax, x_values, y_values, DISPLAY_NAMES[version],
                                    VERSION_COLORS[version])
                        has_data = True

        if has_data:
            ax.set_xlabel('Number of Genomes', fontweight='bold')
            ax.set_ylabel('Context Switches (×10³)', fontweight='bold')
            ax.set_title('(C) I/O: Context Switching', fontweight='bold', pad=10)
            ax.set_xscale('log', base=2)
            ax.set_yscale('log', base=2)
            setup_log_axis_labels(ax, 'y')
            ax.grid(True, alpha=0.3)
            ax.set_xticks([6400, 12800, 25600, 51200, 102400, 204800])
            ax.set_xticklabels(['6400', '12800', '25600', '51200', '102400', '204800'])
        else:
            ax.text(0.5, 0.5, 'No context switch data', ha='center', va='center',
                   fontsize=12, fontweight='bold', transform=ax.transAxes)
            ax.set_title('(C) I/O: Context Switching', fontweight='bold', pad=10)

        for label in ax.get_xticklabels() + ax.get_yticklabels():
            label.set_fontweight('bold')

    # Panel D:
    ax = fig.add_subplot(gs[1, 1])
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
            ax.set_title('(D) I/O: System vs User Time', fontweight='bold', pad=10)
            ax.set_xscale('log', base=2)
            ax.grid(True, alpha=0.3)
            ax.set_xticks([6400, 12800, 25600, 51200, 102400, 204800])
            ax.set_xticklabels(['6400', '12800', '25600', '51200', '102400', '204800'])
        else:
            ax.text(0.5, 0.5, 'No timing breakdown data', ha='center', va='center',
                   fontsize=12, fontweight='bold', transform=ax.transAxes)
            ax.set_title('(D) I/O: System vs User Time', fontweight='bold', pad=10)

        for label in ax.get_xticklabels() + ax.get_yticklabels():
            label.set_fontweight('bold')

    # Add unified legend in the third column
    ax_legend = fig.add_subplot(gs[:, 2])
    ax_legend.axis('off')

    legend = ax_legend.legend(legend_lines, legend_labels,
                              loc='center', fontsize=12,
                              title='Version', title_fontsize=13,
                              frameon=True, fancybox=True, shadow=True)

    legend.get_title().set_fontweight('bold')
    plt.setp(legend.get_texts(), fontweight='bold')

    plt.tight_layout()
    fig.savefig(f'{output_dir}/FigureS1.jpg', dpi=300, bbox_inches='tight')
    plt.close()
    print("Created Figure 1: Rearranged panels (2x2 layout)")

def create_figure2(output_dir, test_data):
    """Create Figure 2: 3x2 layout"""
    fig = plt.figure(figsize=(14, 15))

    # Create gridspec with 3x3 layout (3x2 plots + legend column)
    gs = gridspec.GridSpec(3, 3, figure=fig, hspace=0.3, wspace=0.25,
                          width_ratios=[1, 1, 0.15])

    all_versions = ['AHv1', 'AHv2', 'AHv3', 'AHv4', 'AHv5']

    # Store legend elements
    legend_lines = []
    legend_labels = []

    # Row 1, Col 1: Memory vs N-bases
    ax = fig.add_subplot(gs[0, 0])
    df2 = test_data[2]

    for version in all_versions:
        version_df = df2[df2['version'] == version]
        if not version_df.empty and 'max_memory_kb' in version_df.columns:
            version_df = version_df.copy()
            version_df['max_memory_gb'] = version_df['max_memory_kb'] / (1024 * 1024)
            x_values = sorted(version_df['n_bases'].unique())
            y_values = [version_df[version_df['n_bases'] == x]['max_memory_gb'].values for x in x_values]

            _, _, _, line = plot_with_ci(ax, x_values, y_values, DISPLAY_NAMES[version],
                                         VERSION_COLORS[version])
            if len(legend_lines) < 5:
                legend_lines.append(line)
                legend_labels.append(DISPLAY_NAMES[version])

    ax.set_xlabel('Number of N-bases', fontweight='bold')
    ax.set_ylabel('Peak Memory (GB)', fontweight='bold')
    ax.set_title('(A) Memory: Primer Degeneracy', fontweight='bold', pad=10)
    ax.grid(True, alpha=0.3)
    ax.set_xticks([0, 2, 4, 6])

    for label in ax.get_xticklabels() + ax.get_yticklabels():
        label.set_fontweight('bold')

    # Row 1, Col 2: Memory vs Mismatches
    ax = fig.add_subplot(gs[0, 1])
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
    ax.set_title('(B) Memory: Mismatch Tolerance', fontweight='bold', pad=10)
    ax.grid(True, alpha=0.3)
    ax.set_xticks(range(7))

    for label in ax.get_xticklabels() + ax.get_yticklabels():
        label.set_fontweight('bold')

    # Row 2, Col 1: Thread scaling efficiency
    ax = fig.add_subplot(gs[1, 0])
    df4 = test_data[4]

    max_efficiency = 0
    for version in all_versions:
        version_df = df4[df4['version'] == version]
        if not version_df.empty:
            x_values = sorted(version_df['threads'].unique())

            runtimes = []
            for x in x_values:
                mean_runtime = version_df[version_df['threads'] == x]['real_time'].mean()
                runtimes.append(mean_runtime)

            if len(runtimes) > 0 and runtimes[0] > 0:
                speedup = [runtimes[0] / r if r > 0 else 0 for r in runtimes]
                efficiency = [speedup[i] / x_values[i] * 100 for i in range(len(speedup))]

                ax.plot(x_values, efficiency, color=VERSION_COLORS[version],
                       label=DISPLAY_NAMES[version], linewidth=2.5)
                max_efficiency = max(max_efficiency, max(efficiency))

    ax.axhline(y=100, color='black', linestyle='--', alpha=0.5, label='Perfect Scaling')
    ax.set_xlabel('Number of Threads', fontweight='bold')
    ax.set_ylabel('Parallel Efficiency (%)', fontweight='bold')
    ax.set_title('(C) Efficiency: Thread Scaling', fontweight='bold', pad=10)
    ax.set_xscale('log', base=2)
    ax.grid(True, alpha=0.3)
    ax.set_ylim(0, min(150, max_efficiency * 1.1))

    for label in ax.get_xticklabels() + ax.get_yticklabels():
        label.set_fontweight('bold')

    # Row 2, Col 2: CPU utilization by thread count
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
            ax.set_title('(D) I/O: CPU Utilization vs Threads', fontweight='bold', pad=10)
            ax.set_xscale('log', base=2)
            ax.grid(True, alpha=0.3)
        else:
            ax.text(0.5, 0.5, 'No CPU utilization data', ha='center', va='center',
                   fontsize=12, fontweight='bold', transform=ax.transAxes)
            ax.set_title('(D) I/O: CPU Utilization vs Threads', fontweight='bold', pad=10)

        for label in ax.get_xticklabels() + ax.get_yticklabels():
            label.set_fontweight('bold')

    # Row 3, Col 1: File buffering efficiency
    ax = fig.add_subplot(gs[2, 0])
    df5 = test_data[5]

    if not df5.empty:
        buffer_speedup = []

        for version in all_versions:
            cold_times = df5[(df5['version'] == version) & (df5['cache_state'] == 'cold')]['real_time'].values
            hot_times = df5[(df5['version'] == version) & (df5['cache_state'] == 'hot')]['real_time'].values

            if len(cold_times) > 0 and len(hot_times) > 0:
                cold_mean = np.mean(cold_times)
                hot_mean = np.mean(hot_times)
                speedup = cold_mean / hot_mean if hot_mean > 0 else 1

                buffer_speedup.append({
                    'version': version,
                    'speedup': speedup,
                    'cold': cold_mean,
                    'hot': hot_mean
                })

        if buffer_speedup:
            x = np.arange(len(buffer_speedup))
            speedups = [d['speedup'] for d in buffer_speedup]
            speedup_errors = []

            for version in all_versions:
                cold_times = df5[(df5['version'] == version) & (df5['cache_state'] == 'cold')]['real_time'].values
                hot_times = df5[(df5['version'] == version) & (df5['cache_state'] == 'hot')]['real_time'].values

                if len(cold_times) > 0 and len(hot_times) > 0:
                    ratios = []
                    for _ in range(1000):
                        cold_sample = np.random.choice(cold_times, len(cold_times), replace=True)
                        hot_sample = np.random.choice(hot_times, len(hot_times), replace=True)
                        ratios.append(np.mean(cold_sample) / np.mean(hot_sample))
                    ci_lower, ci_upper = np.percentile(ratios, [2.5, 97.5])
                    speedup_errors.append((ci_upper - ci_lower) / 2)
                else:
                    speedup_errors.append(0)

            colors = [VERSION_COLORS[d['version']] for d in buffer_speedup]

            bars = ax.bar(x, speedups, yerr=speedup_errors[:len(speedups)],
                          color=colors, alpha=0.8, capsize=5,
                          error_kw={'linewidth': 1.5, 'elinewidth': 1.5})

            ax.axhline(y=1, color='black', linestyle='--', alpha=0.5, label='No Buffering Benefit')
            ax.set_xlabel('Version', fontweight='bold')
            ax.set_ylabel('Buffering Speedup Factor', fontweight='bold')
            ax.set_title('(E) Efficiency: File Buffering Benefit', fontweight='bold', pad=10)
            ax.set_xticks(x)
            ax.set_xticklabels([DISPLAY_NAMES[d['version']] for d in buffer_speedup], fontweight='bold')
            ax.grid(True, alpha=0.3, axis='y')
            ax.set_ylim(0, max(speedups) * 1.1)

            for label in ax.get_yticklabels():
                label.set_fontweight('bold')

    # Row 3, Col 2: Speedup comparison
    ax = fig.add_subplot(gs[2, 1])
    df1 = test_data[1]

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
    ax.set_title('(F) Efficiency: Speedup Relative to AHv1.1', fontweight='bold', pad=10)
    ax.set_xscale('log', base=2)
    ax.grid(True, alpha=0.3)
    ax.set_xticks([6400, 12800, 25600, 51200, 102400, 204800])
    ax.set_xticklabels(['6400', '12800', '25600', '51200', '102400', '204800'])

    for label in ax.get_xticklabels() + ax.get_yticklabels():
        label.set_fontweight('bold')

    # Add unified legend in the third column
    ax_legend = fig.add_subplot(gs[:, 2])
    ax_legend.axis('off')

    legend = ax_legend.legend(legend_lines, legend_labels,
                              loc='center', fontsize=12,
                              title='Version', title_fontsize=13,
                              frameon=True, fancybox=True, shadow=True)

    legend.get_title().set_fontweight('bold')
    plt.setp(legend.get_texts(), fontweight='bold')

    plt.tight_layout()
    fig.savefig(f'{output_dir}/FigureS2.jpg', dpi=300, bbox_inches='tight')
    plt.close()
    print("Created Figure 2: Rearranged panels (3x2 layout)")

def main():
    if len(sys.argv) < 2:
        print("Usage: python generate_publication_figures_v5.py <results_dir> [output_dir]")
        sys.exit(1)

    results_dir = sys.argv[1]
    output_dir = sys.argv[2] if len(sys.argv) > 2 else os.path.join(results_dir, 'publication_plots_v5')

    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)

    print(f"Loading data from: {results_dir}")
    print(f"Saving figures to: {output_dir}")
    print("-" * 50)

    # Collect data for all tests
    test_data = {}
    for test_num in range(1, 7):  # Tests 1-6
        print(f"Loading Test {test_num} data...")
        test_data[test_num] = collect_test_data(results_dir, test_num)
        print(f"  Found {len(test_data[test_num])} data points")

    print("-" * 50)

    # Generate all figures
    create_main_figure(output_dir, test_data)
    create_figure1(output_dir, test_data)
    create_figure2(output_dir, test_data)

    print("-" * 50)
    print("All figures generated successfully!")
    print(f"\nOutput directory: {output_dir}")
    print("\nGenerated files:")
    print("  - Figure1_runtime_comparison.jpg (original 3x2 layout)")
    print("  - FigureS1.jpg (2x2 layout)")
    print("  - FigureS2.jpg (3x2 layout)")

if __name__ == "__main__":
    main()
