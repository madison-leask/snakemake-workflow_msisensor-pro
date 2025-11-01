"""
Standard preamble for Jupyter notebooks in the MSIsensor-pro workflow.

This module provides common imports, settings, and utility functions
for analysis notebooks to ensure consistency and reduce code duplication.

Usage in notebooks:
    from project_utils.notebookpreamble import *
"""

# Standard library imports
import sys
from pathlib import Path
import warnings

# Data manipulation and analysis
import pandas as pd
import numpy as np

# Visualization
import matplotlib.pyplot as plt
import seaborn as sns

# Configure display and plotting settings
# Pandas display options
pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', 100)
pd.set_option('display.precision', 3)
pd.set_option('display.width', None)

# Matplotlib/Seaborn settings
plt.style.use('seaborn-v0_8-darkgrid')
sns.set_palette("husl")
plt.rcParams['figure.figsize'] = (12, 8)
plt.rcParams['figure.dpi'] = 100
plt.rcParams['font.size'] = 11

# Suppress warnings for cleaner notebook output
warnings.filterwarnings('ignore', category=FutureWarning)
warnings.filterwarnings('ignore', category=UserWarning)

# Enable automatic reloading of modules
try:
    from IPython import get_ipython
    ipython = get_ipython()
    if ipython is not None:
        ipython.magic('load_ext autoreload')
        ipython.magic('autoreload 2')
        ipython.magic('matplotlib inline')
except (ImportError, AttributeError):
    pass


def setup_paths(workflow_dir=None):
    """
    Set up paths for the workflow.
    
    Parameters
    ----------
    workflow_dir : str or Path, optional
        Path to the workflow directory. If None, assumes notebook is in
        workflow/notebooks/ directory.
    
    Returns
    -------
    dict
        Dictionary containing common paths: 'workflow', 'results', 'config', 'resources'
    """
    if workflow_dir is None:
        # Assume notebook is in workflow/notebooks/
        notebook_dir = Path.cwd()
        workflow_dir = notebook_dir.parent
    else:
        workflow_dir = Path(workflow_dir)
    
    paths = {
        'workflow': workflow_dir,
        'results': workflow_dir.parent / 'results',
        'config': workflow_dir.parent / 'config',
        'resources': workflow_dir.parent / 'resources',
        'logs': workflow_dir.parent / 'logs',
    }
    
    # Add to sys.path for imports
    if str(workflow_dir) not in sys.path:
        sys.path.insert(0, str(workflow_dir))
    
    return paths


def load_config(config_path=None):
    """
    Load workflow configuration.
    
    Parameters
    ----------
    config_path : str or Path, optional
        Path to config.yaml file. If None, uses default location.
    
    Returns
    -------
    dict
        Configuration dictionary
    """
    import yaml
    
    if config_path is None:
        paths = setup_paths()
        config_path = paths['config'] / 'config.yaml'
    else:
        config_path = Path(config_path)
    
    with open(config_path, 'r') as f:
        config = yaml.safe_load(f)
    
    return config


def get_workflow_mode(config=None):
    """
    Determine workflow mode (tumor_only or tumor_normal) from config.
    
    The workflow mode is determined by which alias is configured:
    - If 'baseline' alias is present: tumor_only mode (uses panel of normals)
    - If 'matched_normal' alias is present: tumor_normal mode (uses matched normals)
    
    Parameters
    ----------
    config : dict, optional
        Configuration dictionary. If None, loads from default location.
    
    Returns
    -------
    str
        Workflow mode: 'tumor_only' or 'tumor_normal'
    """
    if config is None:
        config = load_config()
    
    # Check for baseline (panel of normals) -> tumor_only mode
    # Otherwise assume tumor_normal mode with matched normals
    return 'tumor_only' if 'baseline' in config.get('aliases', {}) else 'tumor_normal'


def get_genome_version(config=None):
    """
    Get genome version string from config.
    
    Parameters
    ----------
    config : dict, optional
        Configuration dictionary. If None, loads from default location.
    
    Returns
    -------
    str
        Genome version string in format: genome.dna.{species}.{build}.{release}
    """
    if config is None:
        config = load_config()
    
    species = config['ref']['species']
    build = config['ref']['build']
    release = config['ref']['release']
    
    return f"genome.dna.{species}.{build}.{release}"


def load_samples(samples_path=None):
    """
    Load sample sheet.
    
    Parameters
    ----------
    samples_path : str or Path, optional
        Path to samples.tsv file. If None, uses default location.
    
    Returns
    -------
    pd.DataFrame
        Sample sheet as a pandas DataFrame
    """
    if samples_path is None:
        paths = setup_paths()
        config = load_config()
        samples_path = paths['config'] / config.get('sample_sheet', 'samples.tsv')
    else:
        samples_path = Path(samples_path)
    
    samples = pd.read_csv(samples_path, sep='\t', dtype={'sample': str})
    return samples


def load_msi_results(results_file):
    """
    Load MSI results from the merged output file.
    
    Parameters
    ----------
    results_file : str or Path
        Path to the merged MSI results TSV file
    
    Returns
    -------
    pd.DataFrame
        MSI results as a pandas DataFrame
    """
    results_file = Path(results_file)
    
    if not results_file.exists():
        raise FileNotFoundError(f"Results file not found: {results_file}")
    
    df = pd.read_csv(results_file, sep='\t')
    return df


def classify_msi_status(msi_score, threshold_high=0.2, threshold_low=0.0):
    """
    Classify MSI status based on MSI score.
    
    Parameters
    ----------
    msi_score : float or pd.Series
        MSI score(s) to classify
    threshold_high : float, default=0.2
        Threshold above which samples are classified as MSI-High
    threshold_low : float, default=0.0
        Threshold below which samples are classified as MSS (microsatellite stable)
    
    Returns
    -------
    str or pd.Series
        MSI status: 'MSI-High', 'MSI-Low', or 'MSS'
    """
    if isinstance(msi_score, (pd.Series, np.ndarray)):
        status = pd.Series(['MSS'] * len(msi_score), index=getattr(msi_score, 'index', None))
        status[msi_score > threshold_low] = 'MSI-Low'
        status[msi_score > threshold_high] = 'MSI-High'
        return status
    else:
        if msi_score > threshold_high:
            return 'MSI-High'
        elif msi_score > threshold_low:
            return 'MSI-Low'
        else:
            return 'MSS'


def plot_msi_distribution(msi_data, score_col='msi_score', group_col=None, 
                          threshold_high=0.2, threshold_low=0.0,
                          title='MSI Score Distribution', figsize=(12, 6)):
    """
    Plot distribution of MSI scores.
    
    Parameters
    ----------
    msi_data : pd.DataFrame
        DataFrame containing MSI scores
    score_col : str, default='msi_score'
        Name of column containing MSI scores
    group_col : str, optional
        Column to group by for separate distributions
    threshold_high : float, default=0.2
        MSI-High threshold (shown as vertical line)
    threshold_low : float, default=0.0
        MSI-Low threshold (shown as vertical line)
    title : str
        Plot title
    figsize : tuple
        Figure size
    
    Returns
    -------
    matplotlib.figure.Figure
        The figure object
    """
    fig, axes = plt.subplots(1, 2, figsize=figsize)
    
    # Histogram
    ax1 = axes[0]
    if group_col:
        for group in msi_data[group_col].unique():
            subset = msi_data[msi_data[group_col] == group]
            ax1.hist(subset[score_col], bins=30, alpha=0.6, label=group)
        ax1.legend()
    else:
        ax1.hist(msi_data[score_col], bins=30, alpha=0.7, color='steelblue')
    
    ax1.axvline(threshold_high, color='red', linestyle='--', 
                label=f'MSI-High threshold ({threshold_high})')
    if threshold_low > 0:
        ax1.axvline(threshold_low, color='orange', linestyle='--',
                   label=f'MSI-Low threshold ({threshold_low})')
    ax1.set_xlabel('MSI Score')
    ax1.set_ylabel('Frequency')
    ax1.set_title(f'{title} - Histogram')
    ax1.legend()
    
    # Box plot
    ax2 = axes[1]
    if group_col:
        msi_data.boxplot(column=score_col, by=group_col, ax=ax2)
        ax2.set_title('')
        plt.suptitle('')
    else:
        ax2.boxplot(msi_data[score_col])
        ax2.set_xticklabels(['All Samples'])
    
    ax2.axhline(threshold_high, color='red', linestyle='--')
    if threshold_low > 0:
        ax2.axhline(threshold_low, color='orange', linestyle='--')
    ax2.set_ylabel('MSI Score')
    ax2.set_title(f'{title} - Box Plot')
    
    plt.tight_layout()
    return fig


def summary_statistics(msi_data, score_col='msi_score'):
    """
    Calculate summary statistics for MSI scores.
    
    Parameters
    ----------
    msi_data : pd.DataFrame
        DataFrame containing MSI scores
    score_col : str, default='msi_score'
        Name of column containing MSI scores
    
    Returns
    -------
    pd.Series
        Summary statistics
    """
    stats = msi_data[score_col].describe()
    stats['median'] = msi_data[score_col].median()
    return stats


# Export main functions and modules for use with 'from notebookpreamble import *'
__all__ = [
    'pd', 'np', 'plt', 'sns',
    'setup_paths', 'load_config', 'load_samples',
    'get_workflow_mode', 'get_genome_version',
    'load_msi_results', 'classify_msi_status',
    'plot_msi_distribution', 'summary_statistics',
    'Path', 'sys', 'warnings'
]
