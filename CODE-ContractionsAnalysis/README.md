# CODE-ContractionsAnalysis

This repository contains the Python script (as a Jupyter notebook) for analyzing yolk contractions in fish embryos. Output of the script includes plot of normalized timeseries, plot of detrended timeseries, plot of instantaneous period and wavelet power, and plot of phase as heatmaps. 

## Sample data
A .csv file with sample data is also included. Sample data contains timeseries of mean pixel value of embryos of goldfish, carp, and zebrafish over time. Mean pixel value of fish embryos is used as proxy of projected area of the yolk, and monitoring it over time allows quantification of yolk contractions.

## Dependencies
This repository also contains scripts developed by Gregor Mönke (https://github.com/tensionhead) for wavelet analysis: `lib_psm.py` and `wavelet_analysis.py`.
All files in this repository must be placed in one directory for the script to run smoothly.

## License
Shield: [![CC BY-NC-SA 4.0][cc-by-nc-sa-shield]][cc-by-nc-sa]

This work is licensed under a
[Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License][cc-by-nc-sa].

[![CC BY-NC-SA 4.0][cc-by-nc-sa-image]][cc-by-nc-sa]

[cc-by-nc-sa]: http://creativecommons.org/licenses/by-nc-sa/4.0/
[cc-by-nc-sa-image]: https://licensebuttons.net/l/by-nc-sa/4.0/88x31.png
[cc-by-nc-sa-shield]: https://img.shields.io/badge/License-CC%20BY--NC--SA%204.0-lightgrey.svg

