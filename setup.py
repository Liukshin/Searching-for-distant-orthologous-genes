from setuptools import setup, find_packages
import os

# Read README file
def read_readme():
    try:
        with open("README.md", "r", encoding="utf-8") as fh:
            return fh.read()
    except FileNotFoundError:
        return "PIHMMI - Protein sequence analysis library"

# Read requirements
def read_requirements():
    try:
        with open("requirements.txt", "r", encoding="utf-8") as fh:
            return [line.strip() for line in fh if line.strip() and not line.startswith("#")]
    except FileNotFoundError:
        return [
            "biopython>=1.79",
            "numpy>=1.21.0",
            "pandas>=1.3.0",
            "scipy>=1.7.0",
            "matplotlib>=3.4.0",
            "requests>=2.25.0",
            "scikit-learn>=1.0.0",
            "phytreeviz>=0.1.0",
        ]

setup(
    name="pihmmi",
    version="1.6.0",
    author="Maksim Liukshin",
    author_email="247034@vut.cz",
    description="Library for protein sequence analysis and distant ortholog search",
    long_description=read_readme(),
    long_description_content_type="text/markdown",
    url="https://github.com/Liukshin/Searching-for-distant-orthologous-genes",
    packages=find_packages(),
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    python_requires=">=3.9",
    install_requires=read_requirements(),
    extras_require={
        "hmm": ["pyhmmer>=0.10.0"],  # Optional dependency for HMM
        "orthodb": [],  # Placeholder for OrthoDB (manual installation required)
        "all": ["pyhmmer>=0.10.0"],
        "dev": [
            "pytest>=6.0",
            "pytest-cov",
            "black",
            "flake8",
            "sphinx",
        ],
    },
    include_package_data=True,
    package_data={
        "pihmmi": ["data/*.fasta"],
    },
    entry_points={
        "console_scripts": [
            "pihmmi=pihmmi.pipelines:main_cli",
        ],
    },
)