# Analysis Notebooks

This directory contains Jupyter notebooks for analyzing MSIsensor-pro workflow results.

## Overview

The notebooks in this directory provide standardized analysis and visualization of microsatellite instability (MSI) detection results and related biomarkers.

## Notebooks

### 1. `msi_analysis.ipynb`
**Purpose**: Core MSI analysis and visualization

**Features**:
- Load and explore MSI detection results from MSIsensor-pro
- Classify samples by MSI status (MSI-High, MSI-Low, MSS)
- Generate summary statistics
- Create publication-quality visualizations
- Export annotated results

**When to use**: After running the MSIsensor-pro workflow to analyze and visualize your MSI results.

### 2. `tmb_msi_combined_analysis.ipynb`
**Purpose**: Combined MSI and TMB (Tumor Mutational Burden) analysis

**Features**:
- Compare MSI and TMB biomarkers
- Analyze correlation between MSI and TMB
- Predict immunotherapy response based on both biomarkers
- Comprehensive visualizations

**Note**: This notebook currently uses simulated TMB data for demonstration. To use real TMB data, you need to add variant calling and TMB calculation to the workflow.

**When to use**: When you have both MSI and TMB data available and want to analyze them together.

## Getting Started

### Prerequisites

1. Complete the MSIsensor-pro workflow analysis
2. Have Python 3.7+ with Jupyter installed
3. Required packages (automatically loaded via the preamble):
   - pandas
   - numpy
   - matplotlib
   - seaborn
   - PyYAML
   - scipy (for TMB notebook)

### Installation

If using conda (recommended):
```bash
# Create environment from workflow
conda env create -f ../envs/tidyverse.yaml
conda activate <env_name>

# Or install additional packages
conda install jupyter pandas numpy matplotlib seaborn pyyaml scipy
```

### Running the Notebooks

1. Navigate to the notebooks directory:
```bash
cd workflow/notebooks
```

2. Start Jupyter:
```bash
jupyter notebook
```

3. Open the desired notebook in your browser

4. Run cells sequentially (Shift+Enter) or use "Run All" from the Cell menu

## Using the Notebook Preamble

All notebooks use a standardized preamble from `project_utils.notebookpreamble.py` that provides:

- **Common imports**: pandas, numpy, matplotlib, seaborn
- **Configured settings**: Display options, plotting styles
- **Utility functions**: 
  - `setup_paths()` - Configure workflow paths
  - `load_config()` - Load workflow configuration
  - `load_samples()` - Load sample sheet
  - `load_msi_results()` - Load MSI results
  - `classify_msi_status()` - Classify MSI status
  - `plot_msi_distribution()` - Visualize MSI distributions
  - `summary_statistics()` - Calculate summary stats

### Example Usage

```python
# At the start of any notebook
from project_utils.notebookpreamble import *

# Setup paths and load data
paths = setup_paths()
config = load_config()
samples = load_samples()

# Load and analyze MSI results
msi_results = load_msi_results(paths['results'] / 'tumor_only.genome.dna.homo_sapiens.GRCh38.114.all_samples.tsv')
msi_results['status'] = classify_msi_status(msi_results['msi_score'])

# Visualize
plot_msi_distribution(msi_results)
```

## Customization

### Modifying the Preamble

To add custom functions or change default settings, edit:
```
workflow/project_utils/notebookpreamble.py
```

Changes will be automatically available in all notebooks that import the preamble.

### Creating New Notebooks

When creating new analysis notebooks:

1. Start with the preamble import:
```python
from project_utils.notebookpreamble import *
```

2. Follow the established structure:
   - Introduction/overview
   - Setup paths and configuration
   - Load data
   - Analysis
   - Visualizations
   - Export results
   - Summary

3. Use descriptive markdown cells to explain your analysis

## Output Files

Notebooks may generate output files in the `results/` directory:
- `*.annotated_results.tsv` - MSI results with added classifications
- `*.msi_tmb_combined.tsv` - Combined MSI/TMB analysis (if applicable)
- Figures can be saved using `plt.savefig()` if needed

## Best Practices

1. **Always run the full workflow before analysis**: Ensure all input files are up-to-date
2. **Document changes**: Add markdown cells explaining any modifications
3. **Version control**: Commit notebooks with cleared outputs using:
   ```bash
   jupyter nbconvert --clear-output --inplace *.ipynb
   ```
4. **Reproducibility**: Document package versions and any manual data manipulations
5. **Reuse the preamble**: Don't duplicate code; add new functions to the preamble

## Troubleshooting

### Import Errors
If you get `ModuleNotFoundError: No module named 'project_utils'`:
- Ensure you're running the notebook from the `workflow/notebooks/` directory
- The preamble automatically adds the workflow directory to the Python path

### File Not Found Errors
If data files are not found:
- Check that the workflow has completed successfully
- Verify paths using `paths = setup_paths()` and inspecting the paths dictionary
- Ensure you're using the correct workflow mode and genome version

### Missing Dependencies
Install missing packages:
```bash
conda install <package_name>
# or
pip install <package_name>
```

## Contributing

When adding new notebooks:
1. Use the standardized preamble
2. Follow the existing notebook structure
3. Add comprehensive documentation
4. Update this README
5. Test with the example data in `.test/`

## References

- [MSIsensor-pro documentation](https://github.com/xjtu-omics/msisensor-pro)
- [Jupyter notebook best practices](https://jupyter-notebook.readthedocs.io/)
- [Snakemake workflow best practices](https://snakemake.readthedocs.io/en/stable/snakefiles/best_practices.html)
