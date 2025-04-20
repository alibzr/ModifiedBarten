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

figure; % Create a new figure
subplot(1, 2, 1); % Create the first subplot in a 1x2 grid
loglog(s_frequency_gabor, sensitivity_gabor, 'LineWidth', 2); % Reduced line width for better visibility side by side
xlabel('Spatial Frequency (cpd)', 'FontSize', 12); % Increased font size
ylabel('Sensitivity', 'FontSize', 12); % Increased font size
ylim([1 500])
title('Spatial CSF for Gabor Stimulus', 'FontSize', 14); % Increased font size
set(gca, 'FontSize', 10); % Increase tick mark font size

subplot(1, 2, 2); % Create the second subplot in a 1x2 grid
loglog(t_frequency_disc, sensitivity_disc, 'LineWidth', 2); % Reduced line width
xlabel('Temporal Frequency (Hz)', 'FontSize', 14); % Increased font size
ylabel('Sensitivity', 'FontSize', 14); % Increased font size
ylim([1 500])
title('Temporal CSF for Disc Stimulus', 'FontSize', 14); % Increased font size
set(gca, 'FontSize', 10); % Increase tick mark font size