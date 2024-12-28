# Description
This repository hosts the companion jupyter notebook for the publication *External Acidity as Performance Descriptor in Polyolefin Cracking using Zeolite-Based Materials* (preprint posted to ChemRxiv under doi.org/10.26434/chemrxiv-2024-4fq4v)
It provides an executable version of the manuscript generating all figures and analyses from raw experimental data.
The raw data is hosted on the Open Science Foundation repository under DOI: https://doi.org/10.17605/OSF.IO/PFXH6

Parts of the Python code provided here will be made availible as seperate packages down in the future.

# Getting started
## In Google Colab
The easiest way to run the notebook is in Google Colab.
Click the link below and run the notebook.
<a href="https://colab.research.google.com/github/MyonicS/External_Acidity_Paper/blob/main/Manuscript_Acidity.ipynb" target="_parent"><img src="https://colab.research.google.com/assets/colab-badge.svg" alt="Open In Colab"/></a>
At release the notebook is running python 3.10. 
Cells at the top of the notebook clone the repository, install packages and download the experimental data to the Colab runtime. Note that these steps can take a couple minutes.  

## Running the notebook locally
The notebook is optimized to run in Google Colab. When running the notebook locally delete or comment out all code cells before importing of modules. 
1. Clone the repository:
   ```sh
   git clone https://github.com/MyonicS/External_Acidity_Paper
   ```
2. (optional but recommended: Create a new python 3.10 environment)
3. Install the [spectrochempy library](https://www.spectrochempy.fr/latest/gettingstarted/install/install_win.html) either by following the instruction of the link or by running the command below:
   ```sh
   pip install spectrochempy
   ```
4. Install the remaining requirements by navigating to the repository and running.
   ```sh
   pip install -r requirements.txt
   ```
1. Download the experimental data from the OSF repository.
   You can do this manually by downloading the ZIP file from [here](https://osf.io/pfxh6/) and unzipping the folder into the repository, or by using the [datahugger](https://github.com/J535D165/datahugger) library
6. Run the notebook in your IDE of choice. We recommend VS code as it allows for interactive plots in-line, which can be enabled using 
   ```python
   %matplotlib widget
   ```
	at the start of a cell.
	
# Citing

When utilizing code or data from this study in an academic publication please cite the following preprint
> Rejman, S. _et al._ External Acidity as Performance Descriptor in Polyolefin Cracking using Zeolite-Based Materials. _ChemRxiv_ (2024). doi:10.26434/chemrxiv-2024-4fq4v

Alternatively, the data itself can be cited as 
> Rejman, S. et al. Experimental data supporting: ‘External Acidity as Performance Descriptor in Polyolefin Cracking using Zeolite-Based Materials’. OSF http://dx.doi.org/10.17605/OSF.IO/PFXH6TH (2024)

and the code can be cited under
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.14563455.svg)](https://doi.org/10.5281/zenodo.14563455)

# Bugs and Comments
Feels free to submit bug reports, comments or questions as issues.

