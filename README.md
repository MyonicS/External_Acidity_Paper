# Description
This repository hosts the companion jupyter notebook for the publication *External Acidity as Performance Descriptor in Polyolefin Cracking using Zeolite-Based Materials* (preprint posted to ChemRxiv under DOI.)
It provides an executable version of the manuscript generating all figures and analyses from raw experimental data.
The raw data is hosted on the Open Science Foundation repository under DOI:

# Getting started
## In Google Collab
The easiest way to run the notebook is in Google Collab.
Click the link below and run the notebook.
<a href="https://colab.research.google.com/github/MyonicS/paper_test/blob/main/example_notebook.ipynb" target="_parent"><img src="https://colab.research.google.com/assets/colab-badge.svg" alt="Open In Colab"/></a>
At release the notebook is running python 3.10. 
Cells at the top of the notebook clone the repository, install packages and download the experimental data to the Collab runtime. Note that these steps can take a couple minutes.  

## Running the notebook locally
The notebook is optimized to run in Google Collab. When running the notebook locally delete or comment out all code cells before importing of modules. 
1. Clone the repository:
   ```sh
   git clone LINK
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
   You can do this manually by downloading the ZIP file from [here](https://osf.io/dashboard) and unzipping the folder into the repository, or by using the [datahugger](https://github.com/J535D165/datahugger) library
6. Run the notebook in your IDE of choice. We recommend VS code as it allows for interactive plots in-line, which can be enabled using 
   ```python
   %matplotlib widget
   ```
	at the start of a cell.
	
# Citing

When utilizing code or data from this study in an academic publication please cite the following preprint:
> Sebastian Rejman, *ChemRxiv*

# Bugs and Comments
Feels free to submit bug reports, comments or questions as issues.

