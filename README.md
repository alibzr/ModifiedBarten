# Spatiotemporal Barten CSF (High Frequency Extension) - MATLAB and Python Implementations

This repository contains implementations of the Barten Contrast Sensitivity Function (CSF) model, modified to better predict contrast sensitivity at high temporal frequencies. The model is based on the work described in:

Bozorgian, A., Ashraf, M., & Mantiuk, R. (2024). Spatiotemporal contrast sensitivity functions: predictions for the critical flicker frequency. *Electronic Imaging, 36*, 1-8.

The repository includes the original MATLAB code and a Python translation.

## About the Model

The `CSF_Barten_HF` model extends Barten's CSF to improve its sensitivity predictions both in spatial and temporal domain, with a focus on critical flicker frequency (CFF).

The model incorporates several key factors influencing visual sensitivity:

* Spatial and temporal frequency 
* Eccentricity
* Luminance adaptation level
* Variations in pupil size
* Distribution of retinal ganglion cell density
* Photon noise and internal neural noise

## Implementations

This repository provides two implementations:

1.  **MATLAB:** The original implementation as a MATLAB class file (`CSF_Barten_HF.m`).
2.  **Python:** A translation of the MATLAB code into a Python class (`csf_barten_hf.py`), using NumPy for numerical operations, aiming for similar functionality and maintaining documentation.

## Installation

### MATLAB

1.  Download or clone this repository.
2.  Add the directory containing `CSF_Barten_HF.m` to your MATLAB path.

No additional toolboxes or external dependencies beyond a standard MATLAB installation are required.

### Python

1.  Ensure you have Python (3.6+) installed.
2.  Install the NumPy library:
    ```bash
    pip install numpy
    ```
3.  Download or clone this repository.
4.  Place the `csf_barten_hf.py` file in your project directory or a location included in your Python path where you can import it.

## Usage

Both implementations provide a class (`CSF_Barten_HF` in MATLAB, `CSF_Barten_HF` in Python) with methods to calculate contrast sensitivity for different stimulus types.

### Basic Usage (MATLAB)

```matlab
clear; close; clc;

% Define spatial frequencies
s_frequency = logspace(log10(0.1), log10(100), 500);

% Create a CSF model object
csf_model = CSF_Barten_HF();

% --- Example for Gabor stimuli ---
% Define parameters in a structure
s_frequency_gabor = logspace(log10(0.1), log10(100), 500);
csf_pars_gabor = struct('s_frequency', s_frequency_gabor, 't_frequency', 0, 'eccentricity', 0, ...
                      'luminance', 100, 'ge_sigma', 1, 'vis_field', 0, 'stimulus_type', 'gabor');
% Calculate sensitivity for Gabor
sensitivity_gabor = csf_model.sensitivity(csf_pars_gabor);

% --- Example for Disc stimuli ---
% Define parameters for Disc
t_frequency_disc = logspace(log10(0.1), log10(100), 500);
csf_pars_disc = struct('t_frequency', t_frequency_disc, 'eccentricity', 0, ...
                      'luminance', 100, 'ge_sigma', 1, 'vis_field', 0, 'stimulus_type', 'disc');
% Calculate sensitivity for Disc (using sensitivity_edge)
sensitivity_disc = csf_model.sensitivity_edge(csf_pars_disc);

figure; 
subplot(1, 2, 1);
loglog(s_frequency_gabor, sensitivity_gabor, 'LineWidth', 2);
xlabel('Spatial Frequency (cpd)', 'FontSize', 14); 
ylabel('Sensitivity', 'FontSize', 14);
ylim([1 500])
title('Spatial CSF for Gabor Stimulus', 'FontSize', 14); 
set(gca, 'FontSize', 10);

subplot(1, 2, 2);
loglog(t_frequency_disc, sensitivity_disc, 'LineWidth', 2);
xlabel('Temporal Frequency (Hz)', 'FontSize', 14);
ylabel('Sensitivity', 'FontSize', 14);
ylim([1 500])
title('Temporal CSF for Disc Stimulus', 'FontSize', 14);
set(gca, 'FontSize', 10);