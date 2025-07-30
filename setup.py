from setuptools import setup, find_packages
import os


def read_readme():
    with open("README.md", "r", encoding="utf-8") as fh:
        return fh.read()

def read_requirements():
    with open("requirements.txt", "r", encoding="utf-8") as fh:
        return [line.strip() for line in fh if line.strip() and not line.startswith("#")]


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
        "hmm": ["pyhmmer>=0.10.0"],
        "orthodb": [],
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

