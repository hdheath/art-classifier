#!/usr/bin/env python3
"""Setup script for flAIr."""

from setuptools import setup, find_packages
from pathlib import Path

# Read README
this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text()

setup(
    name="flAIr",
    version="0.1.0",
    description="Functional Long-read Analyzer with Intelligent Reasoning",
    long_description=long_description,
    long_description_content_type="text/markdown",
    author="Harrison Heath",
    author_email="hdheath@ucsc.edu",
    url="https://github.com/yourusername/flAIr",
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    python_requires=">=3.9",
    install_requires=[
        "numpy>=1.24.0",
        "pandas>=2.0.0",
        "scipy>=1.10.0",
        "pyarrow>=12.0.0",
        "torch>=2.0.0",
        "pysam>=0.21.0",
        "pyranges>=0.0.120",
        "matplotlib>=3.7.0",
        "seaborn>=0.12.0",
        "tqdm>=4.65.0",
        "pyyaml>=6.0",
    ],
    extras_require={
        "dev": [
            "pytest>=7.3.0",
            "pytest-cov>=4.1.0",
            "black>=23.3.0",
            "flake8>=6.0.0",
            "mypy>=1.3.0",
        ],
        "docs": [
            "sphinx>=6.2.0",
            "sphinx-rtd-theme>=1.2.0",
        ],
    },
    entry_points={
        "console_scripts": [
            "flair-standardize=extract.standardize_priors:main",
            "flair-extract=extract.bam_to_reads:main",
            "flair-infer=infer.infer_genome:main",
            "flair-eval=eval.compare_baselines:main",
        ],
    },
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
    ],
    keywords="transcriptomics long-read sequencing bioinformatics",
    project_urls={
        "Bug Reports": "https://github.com/yourusername/flAIr/issues",
        "Source": "https://github.com/yourusername/flAIr",
    },
)
