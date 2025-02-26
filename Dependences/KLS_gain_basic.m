function [conversion, offset] = KLS_gain_basic(gain)
% This function describes the fit of conversion factor (y) v. gain setting
    % (x) with an exponential function. When gain is 0 the camera operates
    % as a CCD instead of EMCCD and thus needs a seperate conversion and
    % offset.
% When generating images for the above data fit make sure the intensities 
    % are sub-saturating for valid gain conversion measurement.
%
% Input = gain setting
% Output = [conversion, offset]
% conversion in units e^-/photon
% offset in units ADU
% PhotonData = (RawData - offset) .* conversion;

    if gain == 0
        conversion = 4.2277; % e^-/photon
        offset = 491; % ADU
    else
        conversion = (2.346)*gain^(-0.9559);
        offset = (-0.02152)*gain+(469.6);
    end
end