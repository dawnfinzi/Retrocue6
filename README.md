# Retrocue6

This directory contains the analysis scripts accompanying the manuscript:
Finzi, Wagner, Itthipuripat & Aron, 2017, Stopping cognition: unexpected events recruit prefrontal inhibitory control and interrupt working memory (submitted)

Analysis and scripts adapted from:
Wagner, J., Wessel, J., Gharemahni, A. & Aron, A.R. 2017, Establishing a right frontal beta signature for stopping action in scalp EEG: implications for testing inhibitory control in other task contexts.

For more details on data and the paradigm please refer to the manuscript (Finzi, Wagner, Itthipuripat & Aron, submitted) when available or contact me at: rebecca.d.finzi@gmail.com. Additional information and the behavioral/EEG data is available [on OSF](https://osf.io/h96ny/).


## Scripts/analysis stream
**1. preproc.m**
	- preprocessing of stop data
	- preprocessing of WM data
	- merging of the two datasets
**2. add_error_codes.m**
	- add response error behavioral data to EEG data
**3. run_ica.m**
	- run ica on the merged datasets
**4. post_ICA.m**
	- apply weights to EEG data
	- compute dipoles
**5. brain_comp_select.m**
	- select the ICs that appear to reflect brain activity
**6. create_study.m**
	- select only the brain components (determined in #5)
	- epoch datasets
	- create eeglab STUDY
	- precompute component measures for clustering
**7. clustering.m**
	- cluster components
	- save right frontal beta cluster

### Right frontal beta cluster
**8. computeERSP.m**
	- compute event related spectral perturbations for right frontal beta clustering
**9. AV_ERSP.m**
	- average cluster ERSP images for each condition, compute significance and plot differences
**10. computeERSP_st.m**
	- compute single trial ERSPs for each subject in right frontal beta cluster
**11. regression.m (calls regressERSP_Z.m)**
	- run regression for each subject testing effect of power in time frequency space on behavioral response error
	- compute significance on group level for standard and novel conditions using standardized beta coefficients
	- plot significant results for each condition
	- for novel condition, pool across subjects and plot power (at peak)  response error for visualization purposes

