# Variant Calling Analysis  

**Author:** Luis Aguilera

[![License](https://img.shields.io/badge/License-BSD_3--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)


# Description:
This code is based on [Circulating tumor DNA sequencing in colorectal cancer patients treated with first-line chemotherapy with anti-EGFR.](https://www.nature.com/articles/s41598-021-95345-4). 

This code aims to identify potential variant allele frequency changes that serve as biomarkers for monitoring treatment response and tumor evolution.

<img src="docs/img.png" alt="VCA pipeline" width="900" />


# Pipeline is located in:

**[Source Code](vca_pipeline.ipynb)** - Jupyter notebook containing the VCA pipeline.


# Installation

To install this repository and all its dependencies, we recommend using [Anaconda](https://www.anaconda.com).

* Clone the repository:
```sh
git clone https://github.com/luisub/circulatory_tumor_DNA_analysis.git
```

* Create a virtual environment from the `environment.yml` file and activate it:
```sh
conda env create -f environment.yml
conda activate vca_env
```



