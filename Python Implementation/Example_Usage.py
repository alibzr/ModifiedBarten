from CSF_Barten_HF import CSF_Barten_HF
import numpy as np
import matplotlib.pyplot as plt
from typing import Dict, Union, Literal, TypedDict, Optional

# Example usage
csf_model = CSF_Barten_HF()

# Create spatial frequency vector for Gabor patches
s_frequency_gabor = np.logspace(np.log10(0.1), np.log10(100), 500).reshape(-1, 1)

# Set parameters for Gabor patches
csf_pars_gabor = {
    's_frequency': s_frequency_gabor,
    't_frequency': 0,
    'eccentricity': 0,
    'luminance': 100,
    'ge_sigma': 1,
    'vis_field': 0,
    'stimulus_type': 'gabor'
}

# Calculate sensitivity for Gabor patches
sensitivity_gabor_0hz = csf_model.sensitivity(csf_pars_gabor)

# Create a double logarithmic plot of sensitivity vs spatial frequency
plt.figure(figsize=(10, 6))
plt.loglog(s_frequency_gabor, sensitivity_gabor_0hz, 'k-', linewidth=2, label='Gabor - 0 Hz')
plt.grid(True, which="both", ls="-", alpha=0.6)
plt.xlabel('Spatial Frequency (cpd)', fontsize=12)
plt.ylabel('Contrast Sensitivity', fontsize=12)
plt.ylim(1, 500)
plt.title('Contrast Sensitivity Function (luminance=100 cd/m²)', fontsize=14)

# Add some additional sensitivity curves for different temporal frequencies (Gabor)
t_frequencies_gabor = [2, 8, 16, 32]
for tf in t_frequencies_gabor:
    csf_pars_gabor['t_frequency'] = tf
    sens_gabor = csf_model.sensitivity(csf_pars_gabor)
    plt.loglog(s_frequency_gabor, sens_gabor, label=f'Gabor - {tf} Hz')

plt.legend(title='Temporal Frequency', fontsize=10)
plt.tight_layout()
plt.show()