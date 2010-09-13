
function [IF_interp phi_interp] = interpolate_IP_and_IF(phi_unwrapped,Ts)

IF = [0 phi_unwrapped(2:end)-phi_unwrapped(1:end-1)]/(2*pi*Ts);
IF_good_indexes = find(~isnan(IF));
IF_good_values = IF(~isnan(IF));
IF_interp = interp1(IF_good_indexes,IF_good_values,1:length(IF));

% IF_interp = medfilt1(IF_interp,20);

phi_good_indexes = find(~isnan(phi_unwrapped));
phi_good_values = phi_unwrapped(~isnan(phi_unwrapped));
% phi_interp = interp1(phi_good_indexes,phi_good_values,1:length(phi_unwrapped));

phi_interp = pchip(phi_good_indexes,phi_good_values,1:length(phi_unwrapped));

% phi_interp = medfilt1(phi_interp,20);

% this need to be shape preserving the prevent negative frequencies
% also, it needs to be a spline the prevent discontinuities in the IF
% so it has to be a pchip


% initial and max points = 6 and using previous 3 samples to compute rdm, 5
% imfs, speech signal processing.

% using spline to interp minima and maxima, with 0 and N+1 indexing
% 6.0165e-05 4 interp1
% 8.1090e-05 4 spline
% 5.6189e-05 4 pchip

% using spline to interp minima and maxima, with 1 and N indexing
%  6.3890e-05 4 interp1
%  7.5885e-05 4 spline
%  4.7334e-05 4 pchip

% using pchip to interp minima and maxima, with 1 and N indexing
% 9.6624e-05 4 interp1
% 8.4263e-06 4 spline
% 1.0093e-04 4 pchip

% using pchip to interp minima and maxima, with 0 and N+1 indexing
% 9.6572e-05 4 interp1
% 8.4157e-06 4 spline  <-
% 1.0091e-04 4 pchip
