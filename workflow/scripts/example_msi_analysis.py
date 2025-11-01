#!/usr/bin/env python3
"""
Example script demonstrating the use of project_utils for MSI analysis.

This script shows how to use the utility functions from project_utils.notebookpreamble
in a standalone Python script (not just in notebooks).
"""

import sys
from pathlib import Path

# Add workflow directory to path
workflow_dir = Path(__file__).parent.parent
sys.path.insert(0, str(workflow_dir))

# Import required utilities from the preamble
# Note: For CLI scripts, we import specific functions rather than using 'import *'
# to maintain clarity and avoid importing IPython-specific functionality
from project_utils.notebookpreamble import (
    pd, np, plt,
    setup_paths, load_config, load_samples, load_msi_results,
    classify_msi_status, summary_statistics,
    get_workflow_mode, get_genome_version
)


def main():
    """Main analysis function."""
    print("="*60)
    print("MSI Analysis Example")
    print("="*60)
    
    # Setup paths
    paths = setup_paths(workflow_dir=workflow_dir)
    print("\nWorkflow paths configured:")
    for key, value in paths.items():
        print(f"  {key}: {value}")
    
    # Load configuration
    try:
        config = load_config()
        print(f"\nConfiguration loaded:")
        print(f"  Species: {config.get('ref', {}).get('species')}")
        print(f"  Build: {config.get('ref', {}).get('build')}")
        print(f"  Release: {config.get('ref', {}).get('release')}")
    except Exception as e:
        print(f"\nWarning: Could not load config: {e}")
        return 1
    
    # Load samples
    try:
        samples = load_samples()
        print(f"\nSample information:")
        print(f"  Total samples: {len(samples)}")
        print(f"  Sample types: {samples['alias'].value_counts().to_dict()}")
    except Exception as e:
        print(f"\nWarning: Could not load samples: {e}")
        return 1
    
    # Try to load MSI results
    workflow_mode = get_workflow_mode(config)
    genome_version = get_genome_version(config)
    
    results_file = paths['results'] / f"{workflow_mode}.{genome_version}.all_samples.tsv"
    
    if results_file.exists():
        print(f"\nLoading MSI results from: {results_file}")
        msi_results = load_msi_results(results_file)
        
        # Classify MSI status
        msi_results['msi_status'] = classify_msi_status(msi_results['msi_score'])
        
        # Print summary
        print(f"\nMSI Analysis Results:")
        print(f"  Samples analyzed: {len(msi_results)}")
        print(f"\n  MSI Status Distribution:")
        for status in ['MSI-High', 'MSI-Low', 'MSS']:
            count = (msi_results['msi_status'] == status).sum()
            pct = count / len(msi_results) * 100
            print(f"    {status}: {count} ({pct:.1f}%)")
        
        # Summary statistics
        stats = summary_statistics(msi_results)
        print(f"\n  MSI Score Statistics:")
        print(f"    Mean: {stats['mean']:.3f}")
        print(f"    Median: {stats['median']:.3f}")
        print(f"    Range: {stats['min']:.3f} - {stats['max']:.3f}")
        
        print("\n" + "="*60)
        print("Analysis complete!")
        print("="*60)
        return 0
    else:
        print(f"\nMSI results file not found: {results_file}")
        print("\nPlease run the workflow first to generate results:")
        print("  snakemake --cores 2 --sdm conda")
        return 1


if __name__ == '__main__':
    sys.exit(main())
