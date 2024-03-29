# Neuromelanin
R scripts for the quantification of intracellular and extracellular Neuromelanin (iNM/eNM) [![DOI](https://zenodo.org/badge/742200160.svg)](https://zenodo.org/doi/10.5281/zenodo.10494621)

# Overview
The following repository contains scripts used for the quantification of iNM and eNM
1) NM brightfield images were initially analyzed using the 'TrueAI' feature of [VS200 software](https://www.olympus-lifescience.com/en/solutions-based-systems/vs200/)
2) Data was exported in .xlsx format with tabs "Object Measurements", "Class Measurements" & "ROI Measurements"
3) The "Object Measurements" tab is used for iNM and eNM quantification in which the tab headers should include "Area [µm²]","ROI","Mean (Gray Intensity Value)","Brainbank_ID"
4) Quantification of iNM and eNM can then be performed using the scripts defined in [quant_NM.R](quant_NM.R)
5) Scripts used for regression analysis and plotting of iNM and eNM are defined in [regression_NM.R](regression_NM.R)


# Notes
Determination of modalities; The modalities of the NM Area distribution (log10) were identified by examination of the the modes of kernel density estimator (KDE) as the bandwidth is varied.
The modes for each bandwidth were plotted (Mode Trace; top) and visual examination determined 2 modes (k=2) were prominant (Distribution; bottom). Points indicate where traces are most vertical. 
Scripts used for the determination of modalities are [mode_define.R](mode_define.R)

![image](https://github.com/zchatt/Neuromelanin/assets/30888259/ec751717-97c1-44f9-9882-e35f46d89eb3)

To quantify iNM and eNM granules for each sample in a cohort we perform k-means clustering (k=2) on the log10(Area [µm²]) and set a 95% Confidence Interval (CI) threshold of the normal distribution under each peak.
Scripts used for the quantification of iNM and eNM are [quant_NM.R](quant_NM.R)

![image](https://github.com/zchatt/Neuromelanin/assets/30888259/95cf90f1-5383-4b1c-9be7-065a131d0d96)

# License
This research was funded in whole or in part by Aligning Science Across Parkinson’s (ASAP-020505) through the Michael J. Fox Foundation for Parkinson’s Research (MJFF). For the purpose of open access, the author has applied a CC BY 4.0 public copyright license to all Author Accepted Manuscripts arising from this submission.
