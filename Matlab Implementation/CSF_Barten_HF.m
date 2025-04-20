classdef CSF_Barten_HF
    % CSF_Barten_HF - Spatiotemporal Barten Contrast Sensitivity Function modified for high temporal frequencies
    %
    % This class implements the Barten CSF model with modifications for high temporal frequencies
    % as described in the paper:
    %   "Bozorgian, A., Ashraf, M., & Mantiuk, R. (2024).
    %   Spatiotemporal contrast sensitivity functions: predictions for the critical flicker frequency.
    %   Electronic Imaging, 36, 1-8."
    %
    % The model accounts for:
    % - Spatial and temporal frequency dependencies
    % - Eccentricity effects
    % - Luminance adaptation
    % - Pupil size variations
    % - Retinal ganglion cell density distribution
    % - Photon noise and neural noise
    %
    % Usage example:
    %   % For Gabor stimuli, use sensitivity():
    %   csf_model = CSF_Barten_HF();
    %   s_frequency_gabor = logspace(log10(0.1), log10(100), 500);
    %   csf_pars_gabor = struct('s_frequency', s_frequency_gabor, 't_frequency', 8, 'eccentricity', 0, ...
    %                     'luminance', 100, 'ge_sigma', 1, 'vis_field', 0, 'stimulus_type', 'gabor');
    %   sensitivity_gabor = csf_model.sensitivity(csf_pars_gabor);
    %
    %   % For disc stimuli, use sensitivity_edge():
    %   t_frequency_disc = logspace(log10(0.1), log10(100), 500);
    %   csf_pars_disc = struct('t_frequency', t_frequency_disc, 'eccentricity', 0, ...
    %                     'luminance', 100, 'ge_sigma', 1, 'vis_field', 0, 'stimulus_type', 'disc');
    %   sensitivity_disc = csf_model.sensitivity_edge(csf_pars_disc);

    properties
        par;
        % Constants used in the model calculations
        PHOTON_CONVERSION_FACTOR = 1.240e6;  % Photon conversion factor
        PUPIL_ABERRATION_CONST = 0.08;       % Constant for calculation of sigma based on pupil diameter
    end

    methods
        function obj = CSF_Barten_HF()
            % Constructor for the CSF_Barten_HF model
            %
            % Returns:
            %   obj: Initialized CSF_Barten_HF object with default parameters

            obj.par = obj.get_default_par();
        end

        function S = sensitivity(obj, csf_pars)
            % Calculate the contrast sensitivity for given parameters
            %
            % Parameters:
            %   csf_pars: Structure with the following fields:
            %     s_frequency: Spatial frequency in cycles per degree (cpd)
            %     t_frequency: Temporal frequency in Hz
            %     eccentricity: Retinal eccentricity in degrees
            %     luminance: Average luminance in cd/m²
            %     ge_sigma: Gaussian envelope sigma in degrees
            %     vis_field: Visual field meridian angle in degrees (0=temporal, 99=superior, 180=nasal, 270=inferior)
            %     stimulus_type: Type of stimulus, 'gabor' or 'disc'

            % Validate input parameters
            csf_pars = obj.validate_parameters(csf_pars);

            % Extract parameters
            u = csf_pars.s_frequency;      % Spatial frequency in cpd
            w = csf_pars.t_frequency;      % Temporal frequency in Hz
            e = csf_pars.eccentricity;     % Eccentricity in degrees
            L = csf_pars.luminance;        % Luminance in cd/m²
            X0 = sqrt(pi * csf_pars.ge_sigma.^2); % Angular field size in degrees
            stimulus_type = csf_pars.stimulus_type; % Explicit stimulus type

            % Get model parameters from consolidated structure
            p = obj.par;
            k = p.k;              % Signal to noise ratio
            eta0 = p.eta0;        % Quantum efficiency constant
            sigma0 = p.sigma0;    % Eye MTF constant
            u0 = p.u0;            % Spatial frequency threshold for lateral inhibition
            Phi00 = p.Phi00;      % Neural noise spectral density in fovea
            T = p.T;              % Integration time in seconds
            Xmax0 = p.Xmax0;      % Maximum integration area in fovea
            Nmax = p.Nmax;        % Maximum number of cycles
            tau10 = p.tau10;      % Time constant for cones
            tau20 = p.tau20;      % Time constant for lateral inhibition
            n1 = p.n1;            % Asymptotic slope for cones
            n2 = p.n2;            % Asymptotic slope for lateral inhibition
            k_tem = p.k_tem;      % Temporal meridian constant
            k_nas = p.k_nas;      % Nasal meridian constant

            % Calculate midget ganglion cell density
            Ngm = obj.ganglion_density(csf_pars, "midget");

            % Calculate the maximum integration area as function of eccentricity
            Xmax = Xmax0 * (0.85./(1+(e/4).^2) + 0.15./(1+(e/12).^2)).^-0.5;
            Ymax = Xmax;

            % Calculate quantum efficiency as function of eccentricity
            eta = eta0 * (0.4./(1+(e/7).^2) + 0.48./(1+(e/20).^2) + 0.12);

            % Calculate neural noise spectral density as function of eccentricity and visual field
            alpha = min(1, abs(csf_pars.vis_field-180)/90);
            M_slope = alpha .* k_tem + (1-alpha) .* k_nas;
            Phi0 = Phi00 * (1 + M_slope.*e);

            % Calculate additional spatial parameters
            Y0 = X0;                                    % Tangential field size
            D = sqrt(4/pi * X0.^2);                     % Stimulus diameter (circular assumption)

            % Calculate pupil diameter based on luminance and field size
            d = 5 - 3.*tanh(0.4*log10(L.*(X0.^2/(40^2))));

            % Calculate retinal illuminance in Troland
            E = (pi*(d.^2)/4) .* L .* (1-(d/9.7).^2+(d/12.4).^4);

            % Calculate standard deviation of the line-spread function
            sigmaret = 1./(sqrt(7.2.*sqrt(3)*Ngm/3600));
            sigma00 = sqrt(sigma0^2 - 0.0312^2);
            sigma = sqrt(sigmaret.^2 + sigma00.^2 + (obj.PUPIL_ABERRATION_CONST*d).^2);

            % Calculate optical MTF based on stimulus type
            if strcmp(stimulus_type, "gabor")
                Mopt = exp(-2*pi^2*sigma.^2.*u.^2.*(1/3600));
            elseif strcmp(stimulus_type, "disc")
                 Mopt = 1; % Barten's original model uses Mopt = 1 for disc stimuli
            else
                 error('Invalid stimulus type: %s. Use "gabor" or "disc"', stimulus_type);
            end

            % Calculate receptive field size based on total ganglion cell density
            Ngt = obj.ganglion_density(csf_pars, 'total');
            rf_size = (33163.2/2)./Ngt;

            % Calculate temporal filter time constants
            tau1 = tau10./(1 + 0.55.*log(1 + ((1+D).^0.6).*E./p.t.*rf_size));
            tau2 = tau20./(1 + 0.37.*log(1 + ((1+D/3.2).^5).*E/120));

            % Calculate temporal MTFs
            H1 = 1./(1+(2*pi.*w.*tau1).^2).^(n1/2);  % Cone responses
            H2 = 1./(1+(2*pi.*w.*tau2).^2).^(n2/2);  % Lateral inhibition

            % Calculate lateral inhibition MTF
            F = 1-(1-exp(-(u./u0).^2)).^0.5;  % Lowpass filter
            Mlat = H1.*(1-H2.*F);

            % Calculate photon noise spectral density
            PhiPhoton = 1./(eta.*obj.PHOTON_CONVERSION_FACTOR.*E);

            % Calculate spatial integration parameters
            X = (X0.^-2 + Xmax.^-2 + (((0.5.*X0).^2+4*e.^2)./((0.5.*X0).^2)+e.^2).*(u.^2)./(Nmax.^2)).^-0.5;
            Y = (Y0.^-2 + Ymax.^-2 + (((0.5.*X0).^2)./((0.5.*X0).^2)+e.^2).*(u.^2)./(Nmax.^2)).^-0.5;

            % Calculate contrast sensitivity
            S = sqrt(2) .* (Mopt./(2*k)) .* sqrt((X.*Y.*T)./(PhiPhoton+Phi0./(Mlat.^2)));
        end

        function S = sensitivity_edge(obj, csf_pars)
            % Calculate contrast sensitivity for edge stimuli (discs)
            %
            % Uses Barten's recommendation for handling discs by using
            % fundamental frequency and third harmonic based on disc size
            %
            % Parameters:
            %   csf_pars: Structure with stimulus parameters including 'ge_sigma'
            %
            % Returns:
            %   S: Contrast sensitivity for edge stimulus

            % Ensure stimulus_type is set to 'disc' for this method
            csf_pars.stimulus_type = 'disc';

            % Calculate harmonic frequencies based on disc size
            % Note: The spatial frequency for the disc is determined here based on ge_sigma
            first_harmonic = 1./(sqrt(pi).*2.*csf_pars.ge_sigma);
            third_harmonic = 3./(sqrt(pi).*2.*csf_pars.ge_sigma);

            % Choose appropriate frequency based on disc size
            small_disc = csf_pars.ge_sigma < 2.5;
            large_disc = csf_pars.ge_sigma >= 2.5;
            
            % Assign the calculated spatial frequency based on disc size
            csf_pars.s_frequency = small_disc.*first_harmonic + large_disc.*third_harmonic;

            % Calculate sensitivity using the main sensitivity method with stimulus_type='disc'
            S = obj.sensitivity(csf_pars);

        end


        function Ng = ganglion_density(obj, csf_pars, cell_type)
            % Calculate ganglion cell density based on Watson's formula
            %
            % Parameters:
            %   csf_pars: Structure with stimulus parameters including 'eccentricity' and 'vis_field'
            %   cell_type: Type of ganglion cells ('midget', 'parasol', or 'total')
            %   stimulus_type: Type of stimulus, 'gabor' or 'disc'
            %
            % Returns:
            %   Ng: Ganglion cell density in cells/deg²

            % Extract parameters
            meridian = csf_pars.vis_field;
            e = csf_pars.eccentricity;

            % Density of total retinal ganglion cells in fovea
            Ng0 = 33163.2;

            % Scale factor for decline in midget fraction with eccentricity
            rm = 41.03;

            % Calculate cell fraction based on cell type
            % Assumed that 50% of a cell type form an independent mosaic
            switch cell_type
                case "midget"
                    f0 = 0.8928/2 .* (1+e./rm).^-1;
                case "parasol"
                    % Parasol/midget ratio (rough estimate from Dacey (1992))
                    f0 = (1-0.8928)/2 .* (1+e./rm).^-1 .* (0.006.*e+0.015);
                case "total"
                    f0 = 1/2;
                otherwise
                    error('Unknown cell type: %s. Use "midget", "parasol", or "total"', cell_type);
            end

            % Convert retinal meridians to visual field
            meridian = rem(meridian+180, 360);
            tem = meridian == 0;   % Temporal
            sup = meridian == 90;  % Superior
            nas = meridian == 180; % Nasal
            inf = meridian == 270; % Inferior

            % Parameters for Watson's formula based on visual meridians
            temPar = tem.*[0.9851 1.058 22.14];
            supPar = sup.*[0.9935 1.035 16.35];
            nasPar = nas.*[0.9729 1.084 7.633];
            infPar = inf.*[0.9960 0.9932 12.13];
            visPar = temPar + supPar + nasPar + infPar;

            % Apply Watson's formula for ganglion density
            a = visPar(:,1);
            r2 = visPar(:,2);
            re = visPar(:,3);
            Ng = Ng0 .* f0 .* (a.*(1+(e./r2)).^-2 + (1-a).*exp(-e./re));

        end
    end

    methods (Access = private)
        function csf_pars = validate_parameters(obj, csf_pars)
            % Validate and complete input parameters
            %
            % Parameters:
            %   csf_pars: Input parameter structure
            %
            % Returns:
            %   csf_pars: Validated parameter structure

            % Check for required parameters, including stimulus_type
            required_params = {'luminance', 'ge_sigma', 's_frequency', 't_frequency', 'eccentricity', 'stimulus_type'};
            csf_pars = obj.test_complete_params(csf_pars, required_params, true);

            % Validate stimulus_type
            if ~ismember(csf_pars.stimulus_type, {'gabor', 'disc'})
                 error('Invalid value for stimulus_type: %s. Use ''gabor'' or ''disc''.', csf_pars.stimulus_type);
            end

            return;
        end

        % Removed determine_stimulus_type method

        function pars = test_complete_params(obj, pars, requires, expand)
            % code adopted from: https://github.com/gfxdisp/castleCSF
            % Modified to include 'stimulus_type' in validation

            if nargin < 3
                expand = false;
            end

            valid_names = { ...
                'luminance', 'lms_bkg', 'lms_delta', ...
                's_frequency', 't_frequency', 'orientation', ...
                'area', 'ge_sigma', 'eccentricity', 'vis_field', ...
                'stimulus_type' % Added stimulus_type
            };

            % Validate field names
            fn = fieldnames(pars);
            cur_par = 1;
            color_dim = [];

            for kk = 1:length(fn)
                name = fn{kk};

                if ~ismember(name, valid_names)
                    error('Parameter structure contains unrecognized field ''%s''', name);
                end

                try
                    param = pars.(name);

                    if ismember(name, {'lms_bkg', 'lms_delta'})
                        p_sz = size(param);
                        if p_sz(end) ~= 3
                            error('The last dimension of ''%s'' must have size 3', name);
                        end

                        if ~isempty(color_dim) && color_dim ~= ndims(param)
                            error('LMS colour must be %d-th dimension of your input parameters', color_dim);
                        end

                        color_dim = ndims(param);
                    end

                    % Exclude stimulus_type from the broadcasting check as it's a string
                    if ~strcmp(name, 'stimulus_type')
                         cur_par = cur_par .* param;  % Used to validate broadcasting
                    end


                catch E
                    if strcmp(E.identifier, 'MATLAB:sizeDimensionsMustMatch')
                        error('Parameter %s cannot be broadcasted', name);
                    else
                        rethrow(E);
                    end
                end
            end

            % Infer color dimension if not found
            if isempty(color_dim)
                % Find the first numeric parameter to infer dimensions if no color param
                 numeric_fields = setdiff(fn, {'stimulus_type'});
                 if ~isempty(numeric_fields)
                    color_dim = ndims(pars.(numeric_fields{1})) + 1;
                 else
                     color_dim = 1; % Default if no numeric parameters
                 end
            end

            % Fill in required parameters
            if ismember('luminance', requires) && ~isfield(pars, 'luminance')
                if ~isfield(pars, 'lms_bkg')
                    error('You need to pass either ''luminance'' or ''lms_bkg'' parameter.');
                end
                pars.luminance = CSF_base.last_dim(pars.lms_bkg, 1) + CSF_base.last_dim(pars.lms_bkg, 2);
            end

            if ismember('lms_bkg', requires) && ~isfield(pars, 'lms_bkg')
                if ~isfield(pars, 'luminance')
                    error('You need to pass either ''luminance'' or ''lms_bkg'' parameter.');
                end
                default_lms_weights = [0.6991 0.3009 0.0198];
                pars.lms_bkg = reshape(default_lms_weights, [ones(1, color_dim-1), 3]) .* pars.luminance;
            end

            if ismember('ge_sigma', requires) && ~isfield(pars, 'ge_sigma')
                if ~isfield(pars, 'area')
                    error('You need to pass either ''ge_sigma'' or ''area'' parameter.');
                end
                pars.ge_sigma = sqrt(pars.area / pi);
            end

            if ismember('area', requires) && ~isfield(pars, 'area')
                if ~isfield(pars, 'ge_sigma')
                    error('You need to pass either ''ge_sigma'' or ''area'' parameter.');
                end
                pars.area = pi * pars.ge_sigma.^2;
            end

             % Check if stimulus_type is required and missing
            if ismember('stimulus_type', requires) && ~isfield(pars, 'stimulus_type')
                 error('The ''stimulus_type'' parameter is required.');
            end


            % Set default parameter values (excluding stimulus_type as it's required)
            def_pars = struct( ...
                'eccentricity', 0, ...
                'vis_field', 180, ...
                'orientation', 0, ...
                't_frequency', 0, ...
                'lms_delta', reshape([0.6855 0.2951 0.0194], [ones(1, color_dim-1), 3]) ...
            );

            def_fields = fieldnames(def_pars);
            for kk = 1:length(def_fields)
                field = def_fields{kk};
                if ~isfield(pars, field)
                    pars.(field) = def_pars.(field);
                end
            end

        end
    end

    methods (Static)
        function p = get_default_par()
            % Returns default model parameters
            %
            % Returns:
            %   p: Structure with default parameter values

            p = struct();
            % Model parameters as originally proposed by Barten (1999)
            p.eta0 = 0.03;        % Quantum efficiency constant
            p.Phi00 = 3e-8;       % Neural noise spectral density in fovea
            p.T = 0.1;            % Integration time in seconds
            p.Xmax0 = 12;         % Maximum integration area in fovea
            p.Nmax = 15;          % Maximum number of cycles
            p.n1 = 7;             % Asymptotic slope for cones
            p.n2 = 4;             % Asymptotic slope for lateral inhibition

            % Model parameters from fitting to experimental data
            p.k = 8.36831;        % A parameter proportional to signal to noise ratio
            p.sigma0 = 0.370199;  % Eye MTF constant
            p.u0 = 2.44734;       % Spatial frequency threshold for lateral inhibition
            p.tau10 = 0.0334009;  % Time constant for cones
            p.tau20 = 0.0477005;  % Time constant for lateral inhibition
            p.k_tem = 0.556637;   % Temporal meridian constant
            p.k_nas = 0.171344;   % Nasal meridian constant
            p.t = 1.32776;        % Time constant scaling factor
        end

        function Y = last_dim(X,d)
            cln(1:ndims(X)) = { ':' };
            cln{end} = d;
            Y = X(cln{:});
        end

    end
end