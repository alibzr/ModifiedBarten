clear; close; clc;

% Define spatial frequencies
s_frequency = logspace(log10(0.1), log10(100), 500);

% Initialize CSF model
csf_model = CSF_Barten_HF();

% Set parameters
csf_pars = struct( ...
    's_frequency', s_frequency, ...
    't_frequency', 0, ...
    'eccentricity', 0, ...
    'luminance', 100, ...
    'ge_sigma', 2, ...
    'vis_field', 0, ...
    'stimulus_type', "gabor" ...
);

% Compute sensitivity
sensitivity = csf_model.sensitivity(csf_pars);

% Plotting
figure;
loglog(s_frequency, sensitivity, '-', 'LineWidth', 5);
grid on;
xlabel('Spatial Frequency (cycles/degree)', 'FontSize', 12);
ylabel('Contrast Sensitivity', 'FontSize', 12);
title('Barten CSF - High Frequency', 'FontSize', 14);
ylim([1 500]);

% Improve ticks and appearance
set(gca, 'FontSize', 14);
legend('CSF Curve', 'Location', 'southwest');
