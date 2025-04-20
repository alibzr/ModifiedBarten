import numpy as np
import math

class CSF_Barten_HF:
    """
    CSF_Barten_HF - Spatiotemporal Barten Contrast Sensitivity Function modified for high temporal frequencies

    This class implements the Barten CSF model with modifications for high temporal frequencies
    as described in the paper:
      "Bozorgian, A., Ashraf, M., & Mantiuk, R. (2024).
      Spatiotemporal contrast sensitivity functions: predictions for the critical flicker frequency.
      Electronic Imaging, 36, 1-8."

    The model accounts for:
    - Spatial and temporal frequency dependencies
    - Eccentricity effects
    - Luminance adaptation
    - Pupil size variations
    - Retinal ganglion cell density distribution
    - Photon noise and neural noise

    Usage example:
      # For Gabor stimuli, use sensitivity():
      csf_model = CSF_Barten_HF()
      s_frequency_gabor = np.logspace(np.log10(0.1), np.log10(100), 500)
      csf_pars_gabor = {'s_frequency': s_frequency_gabor, 't_frequency': 8, 'eccentricity': 0,
                        'luminance': 100, 'ge_sigma': 1, 'vis_field': 0, 'stimulus_type': 'gabor'}
      sensitivity_gabor = csf_model.sensitivity(csf_pars_gabor)

      # For disc stimuli, use sensitivity_edge():
      t_frequency_disc = np.logspace(np.log10(0.1), np.log10(100), 500)
      csf_pars_disc = {'t_frequency': t_frequency_disc, 'eccentricity': 0,
                        'luminance': 100, 'ge_sigma': 1, 'vis_field': 0, 'stimulus_type': 'disc'}
      sensitivity_disc = csf_model.sensitivity_edge(csf_pars_disc)
    """

    # Constants used in the model calculations
    PHOTON_CONVERSION_FACTOR = 1.240e6  # Photon conversion factor
    PUPIL_ABERRATION_CONST = 0.08       # Constant for calculation of sigma based on pupil diameter

    def __init__(self):
        """
        Constructor for the CSF_Barten_HF model

        Initializes the CSF_Barten_HF object with default parameters.
        """
        self.par = self.get_default_par()

    def sensitivity(self, csf_pars):
        """
        Calculate the contrast sensitivity for given parameters

        Parameters:
            csf_pars (dict): Dictionary with the following keys:
              s_frequency (float or np.ndarray): Spatial frequency in cycles per degree (cpd)
              t_frequency (float or np.ndarray): Temporal frequency in Hz
              eccentricity (float or np.ndarray): Retinal eccentricity in degrees
              luminance (float or np.ndarray): Average luminance in cd/m²
              ge_sigma (float or np.ndarray): Gaussian envelope sigma in degrees
              vis_field (float or np.ndarray): Visual field meridian angle in degrees (0=temporal, 99=superior, 180=nasal, 270=inferior)
              stimulus_type (str): Type of stimulus, 'gabor' or 'disc'

        Returns:
            np.ndarray: Contrast sensitivity
        """
        # Validate input parameters
        csf_pars = self._validate_parameters(csf_pars)

        # Extract parameters
        u = csf_pars['s_frequency']      # Spatial frequency in cpd
        w = csf_pars['t_frequency']      # Temporal frequency in Hz
        e = csf_pars['eccentricity']     # Eccentricity in degrees
        L = csf_pars['luminance']        # Luminance in cd/m²
        X0 = np.sqrt(np.pi * csf_pars['ge_sigma']**2) # Angular field size in degrees
        stimulus_type = csf_pars['stimulus_type'] # Explicit stimulus type

        # Get model parameters from consolidated dictionary
        p = self.par
        k = p['k']              # Signal to noise ratio
        eta0 = p['eta0']        # Quantum efficiency constant
        sigma0 = p['sigma0']    # Eye MTF constant
        u0 = p['u0']            # Spatial frequency threshold for lateral inhibition
        Phi00 = p['Phi00']      # Neural noise spectral density in fovea
        T = p['T']              # Integration time in seconds
        Xmax0 = p['Xmax0']      # Maximum integration area in fovea
        Nmax = p['Nmax']        # Maximum number of cycles
        tau10 = p['tau10']      # Time constant for cones
        tau20 = p['tau20']      # Time constant for lateral inhibition
        n1 = p['n1']            # Asymptotic slope for cones
        n2 = p['n2']            # Asymptotic slope for lateral inhibition
        k_tem = p['k_tem']      # Temporal meridian constant
        k_nas = p['k_nas']      # Nasal meridian constant

        # Calculate midget ganglion cell density
        Ngm = self.ganglion_density(csf_pars, "midget")

        # Calculate the maximum integration area as function of eccentricity
        Xmax = Xmax0 * (0.85/(1+(e/4)**2) + 0.15/(1+(e/12)**2))**-0.5
        Ymax = Xmax

        # Calculate quantum efficiency as function of eccentricity
        eta = eta0 * (0.4/(1+(e/7)**2) + 0.48/(1+(e/20)**2) + 0.12)

        # Calculate neural noise spectral density as function of eccentricity and visual field
        alpha = np.minimum(1, np.abs(csf_pars['vis_field']-180)/90)
        M_slope = alpha * k_tem + (1-alpha) * k_nas
        Phi0 = Phi00 * (1 + M_slope*e)

        # Calculate additional spatial parameters
        Y0 = X0                                    # Tangential field size
        D = np.sqrt(4/np.pi * X0**2)                     # Stimulus diameter (circular assumption)

        # Calculate pupil diameter based on luminance and field size
        d = 5 - 3*np.tanh(0.4*np.log10(L * (X0**2 / (40**2))))

        # Calculate retinal illuminance in Troland
        E = (np.pi*(d**2)/4) * L * (1-(d/9.7)**2+(d/12.4)**4)

        # Calculate standard deviation of the line-spread function
        sigmaret = 1/(np.sqrt(7.2*np.sqrt(3)*Ngm/3600))
        sigma00 = np.sqrt(sigma0**2 - 0.0312**2)
        sigma = np.sqrt(sigmaret**2 + sigma00**2 + (self.PUPIL_ABERRATION_CONST*d)**2)

        # Calculate optical MTF based on stimulus type
        if stimulus_type == "gabor":
            Mopt = np.exp(-2*np.pi**2*sigma**2*u**2*(1/3600))
        elif stimulus_type == "disc":
             Mopt = 1 # Barten's original model uses Mopt = 1 for disc stimuli
        else:
             raise ValueError(f"Invalid stimulus type: {stimulus_type}. Use 'gabor' or 'disc'")

        # Calculate receptive field size based on total ganglion cell density
        Ngt = self.ganglion_density(csf_pars, 'total')
        # Ensure no division by zero if Ngt can be zero
        rf_size = (33163.2/2) / np.maximum(Ngt, np.finfo(float).eps) # Add small epsilon to avoid division by zero

        # Calculate temporal filter time constants
        # Ensure positive argument for log
        tau1 = tau10 / (1 + 0.55 * np.log(1 + np.maximum(1e-9, ((1+D)**0.6) * E / p['t'] * rf_size)))
        tau2 = tau20 / (1 + 0.37 * np.log(1 + np.maximum(1e-9, ((1+D/3.2)**5) * E / 120)))

        # Calculate temporal MTFs
        H1 = 1 / (1+(2*np.pi*w*tau1)**2)**(n1/2)  # Cone responses
        H2 = 1 / (1+(2*np.pi*w*tau2)**2)**(n2/2)  # Lateral inhibition

        # Calculate lateral inhibition MTF
        # Ensure positive argument for sqrt
        F = 1 - (1 - np.exp(-(u/u0)**2))**0.5  # Lowpass filter
        Mlat = H1 * (1-H2*F)

        # Calculate photon noise spectral density
        # Ensure no division by zero if E or eta can be zero
        PhiPhoton = 1 / (eta * self.PHOTON_CONVERSION_FACTOR * np.maximum(E, np.finfo(float).eps))


        # Calculate spatial integration parameters
        # Ensure positive arguments for power and sqrt
        X = (X0**-2 + Xmax**-2 + (((0.5*X0)**2+4*e**2)/((0.5*X0)**2)+e**2)*(u**2)/(Nmax**2))**-0.5
        Y = (Y0**-2 + Ymax**-2 + (((0.5*X0)**2)/((0.5*X0)**2)+e**2)*(u**2)/(Nmax**2))**-0.5

        # Ensure Mlat is not zero before division
        Mlat_safe = np.maximum(Mlat, np.finfo(float).eps)

        # Calculate contrast sensitivity
        S = np.sqrt(2) * (Mopt/(2*k)) * np.sqrt((X*Y*T)/(PhiPhoton+Phi0/(Mlat_safe**2)))

        return S

    def sensitivity_edge(self, csf_pars):
        """
        Calculate contrast sensitivity for edge stimuli (discs)

        Uses Barten's recommendation for handling discs by using
        fundamental frequency and third harmonic based on disc size

        Parameters:
            csf_pars (dict): Dictionary with stimulus parameters including 'ge_sigma'

        Returns:
            np.ndarray: Contrast sensitivity for edge stimulus
        """
        # Ensure stimulus_type is set to 'disc' for this method
        csf_pars['stimulus_type'] = 'disc'

        # Calculate harmonic frequencies based on disc size
        # Note: The spatial frequency for the disc is determined here based on ge_sigma
        # Ensure ge_sigma is not zero before division
        ge_sigma_safe = np.maximum(csf_pars['ge_sigma'], np.finfo(float).eps)
        first_harmonic = 1 / (np.sqrt(np.pi) * 2 * ge_sigma_safe)
        third_harmonic = 3 / (np.sqrt(np.pi) * 2 * ge_sigma_safe)

        # Choose appropriate frequency based on disc size
        # Use np.where for element-wise conditional selection
        csf_pars['s_frequency'] = np.where(csf_pars['ge_sigma'] < 2.5, first_harmonic, third_harmonic)

        # Calculate sensitivity using the main sensitivity method with stimulus_type='disc'
        S = self.sensitivity(csf_pars)

        return S

    def ganglion_density(self, csf_pars, cell_type):
        """
        Calculate ganglion cell density based on Watson's formula

        Parameters:
            csf_pars (dict): Dictionary with stimulus parameters including 'eccentricity' and 'vis_field'
            cell_type (str): Type of ganglion cells ('midget', 'parasol', or 'total')

        Returns:
            np.ndarray: Ganglion cell density in cells/deg²
        """
        # Extract parameters
        meridian = csf_pars['vis_field']
        e = csf_pars['eccentricity']

        # Density of total retinal ganglion cells in fovea
        Ng0 = 33163.2

        # Scale factor for decline in midget fraction with eccentricity
        rm = 41.03

        # Calculate cell fraction based on cell type
        # Assumed that 50% of a cell type form an independent mosaic
        if cell_type == "midget":
            f0 = 0.8928/2 * (1+e/rm)**-1
        elif cell_type == "parasol":
            # Parasol/midget ratio (rough estimate from Dacey (1992))
            f0 = (1-0.8928)/2 * (1+e/rm)**-1 * (0.006*e+0.015)
        elif cell_type == "total":
            f0 = 1/2
        else:
            raise ValueError(f"Unknown cell type: {cell_type}. Use 'midget', 'parasol', or 'total'")

        # Convert retinal meridians to visual field
        meridian = np.remainder(meridian+180, 360)

        # Parameters for Watson's formula based on visual meridians
        # Using np.where to handle different meridians
        a = np.where(meridian == 0, 0.9851,
              np.where(meridian == 90, 0.9935,
              np.where(meridian == 180, 0.9729,
              np.where(meridian == 270, 0.9960, 0)))) # Default to 0 if not a recognized meridian

        r2 = np.where(meridian == 0, 1.058,
               np.where(meridian == 90, 1.035,
               np.where(meridian == 180, 1.084,
               np.where(meridian == 270, 0.9932, 0))))

        re = np.where(meridian == 0, 22.14,
               np.where(meridian == 90, 16.35,
               np.where(meridian == 180, 7.633,
               np.where(meridian == 270, 12.13, 0))))

        # Apply Watson's formula for ganglion density
        # Ensure no division by zero if r2 can be zero
        r2_safe = np.maximum(r2, np.finfo(float).eps)
        Ng = Ng0 * f0 * (a * (1 + e / r2_safe)**-2 + (1 - a) * np.exp(-e / np.maximum(re, np.finfo(float).eps)))

        return Ng

    def _validate_parameters(self, csf_pars):
        """
        Validate and complete input parameters (private helper method)

        Parameters:
            csf_pars (dict): Input parameter dictionary

        Returns:
            dict: Validated parameter dictionary
        """
        # Check for required parameters, including stimulus_type
        required_params = ['luminance', 'ge_sigma', 's_frequency', 't_frequency', 'eccentricity', 'stimulus_type']
        csf_pars = self._test_complete_params(csf_pars, required_params)

        # Validate stimulus_type
        if csf_pars['stimulus_type'] not in ['gabor', 'disc']:
             raise ValueError(f"Invalid value for stimulus_type: {csf_pars['stimulus_type']}. Use 'gabor' or 'disc'.")

        return csf_pars

    def _test_complete_params(self, pars, requires):
        """
        Complete parameter dictionary with default values and validate (private helper method)

        Code adopted from: https://github.com/gfxdisp/castleCSF
        Modified to include 'stimulus_type' in validation

        Parameters:
            pars (dict): Input parameter dictionary
            requires (list): List of required parameter names

        Returns:
            dict: Completed and validated parameter dictionary
        """
        valid_names = {
            'luminance', 'lms_bkg', 'lms_delta',
            's_frequency', 't_frequency', 'orientation',
            'area', 'ge_sigma', 'eccentricity', 'vis_field',
            'stimulus_type' # Added stimulus_type
        }

        # Validate field names
        for name in pars.keys():
            if name not in valid_names:
                raise ValueError(f"Parameter dictionary contains unrecognized key '{name}'")

        # Infer color dimension if needed (simplified for this model as LMS conversion is not used)
        # In a full implementation handling color, this logic would be more complex.
        # For this model, we assume parameters are scalar or broadcastable NumPy arrays.
        color_dim = 1 # Assuming no explicit color dimension in input for this model

        # Fill in required parameters
        if 'luminance' in requires and 'luminance' not in pars:
            if 'lms_bkg' not in pars:
                raise ValueError('You need to pass either \'luminance\' or \'lms_bkg\' parameter.')
            # Simplified LMS to luminance conversion for validation purpose if needed
            # This part might need adjustment if lms_bkg is actually used in sensitivity calculation
            pars['luminance'] = np.sum(pars['lms_bkg'], axis=-1) # Example if lms_bkg is [..., 3]

        if 'lms_bkg' in requires and 'lms_bkg' not in pars:
            if 'luminance' not in pars:
                 raise ValueError('You need to pass either \'luminance\' or \'lms_bkg\' parameter.')
            default_lms_weights = np.array([0.6991, 0.3009, 0.0198])
            # Reshape for broadcasting if luminance is an array
            pars['lms_bkg'] = default_lms_weights * np.reshape(pars['luminance'], pars['luminance'].shape + (1,))


        if 'ge_sigma' in requires and 'ge_sigma' not in pars:
            if 'area' not in pars:
                raise ValueError('You need to pass either \'ge_sigma\' or \'area\' parameter.')
            pars['ge_sigma'] = np.sqrt(pars['area'] / np.pi)

        if 'area' in requires and 'area' not in pars:
            if 'ge_sigma' not in pars:
                raise ValueError('You need to pass either \'ge_sigma\' or \'area\' parameter.')
            pars['area'] = np.pi * pars['ge_sigma']**2

        # Check if stimulus_type is required and missing
        if 'stimulus_type' in requires and 'stimulus_type' not in pars:
             raise ValueError('The \'stimulus_type\' parameter is required.')

        # Set default parameter values (excluding stimulus_type as it's required)
        def_pars = {
            'eccentricity': 0,
            'vis_field': 180,
            'orientation': 0,
            't_frequency': 0,
            'lms_delta': np.array([0.6855, 0.2951, 0.0194]) # Assuming no color dim for simplicity
        }

        for field, default_value in def_pars.items():
            if field not in pars:
                pars[field] = default_value

        # Ensure all required parameters are present after defaults
        for req in requires:
            if req not in pars:
                 raise ValueError(f"Required parameter '{req}' is missing.")

        return pars

    @staticmethod
    def get_default_par():
        """
        Returns default model parameters (static method)

        Returns:
            dict: Dictionary with default parameter values
        """
        p = {}
        # Model parameters as originally proposed by Barten (1999)
        p['eta0'] = 0.03        # Quantum efficiency constant
        p['Phi00'] = 3e-8       # Neural noise spectral density in fovea
        p['T'] = 0.1            # Integration time in seconds
        p['Xmax0'] = 12         # Maximum integration area in fovea
        p['Nmax'] = 15          # Maximum number of cycles
        p['n1'] = 7             # Asymptotic slope for cones
        p['n2'] = 4             # Asymptotic slope for lateral inhibition
        # Model parameters from fitting to experimental data
        p['k'] = 8.36831        # A parameter proportional to signal to noise ratio
        p['sigma0'] = 0.370199  # Eye MTF constant
        p['u0'] = 2.44734       # Spatial frequency threshold for lateral inhibition
        p['tau10'] = 0.0334009  # Time constant for cones
        p['tau20'] = 0.0477005  # Time constant for lateral inhibition
        p['k_tem'] = 0.556637   # Temporal meridian constant
        p['k_nas'] = 0.171344   # Nasal meridian constant
        p['t'] = 1.32776        # Time constant scaling factor (Note: 't' is also used for integration time in the paper, renamed to avoid conflict)
        return p

    @staticmethod
    def last_dim(X, d):
        """
        Equivalent of MATLAB's X(:,:,d) for the last dimension
        (static method)

        This function is likely not needed with proper NumPy broadcasting
        and array indexing, but included for completeness based on the
        original MATLAB code's helper function.
        """
        # This implementation assumes d is an index for the last dimension
        # and X is a NumPy array.
        if not isinstance(X, np.ndarray):
            raise TypeError("Input X must be a NumPy array.")
        if d < 0 or d >= X.shape[-1]:
             raise IndexError("Index d out of bounds for the last dimension.")

        # Construct a tuple of slices to select the d-th element of the last dimension
        slices = (slice(None),) * (X.ndim - 1) + (d,)
        return X[slices]