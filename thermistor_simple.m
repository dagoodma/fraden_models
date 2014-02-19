function [ R T ] = thermistor_simple( r0, t0, r1, t1 )
%THERMISTOR_SIMPLE Simple model for an NTC thermistor.
%   Given two measured resistances, r0 and r1, and known temperatures in
%   celcius t0 and t1, this function will output two anonymous functions
%   that model the resitance as a function of temperature and visa versa.
%
if (r0 == 0 || t0 == 0 || r1 == 0 || t1 == 0)
    error('All arguments must be non-zero values.');
end

% Celcius to kelvin conversion functions
celcius_to_kelvin = @(TC) TC + 273.15;
kelvin_to_celcius = @(TK) TK - 273.15;

% Convert calibration points from celcius to kelvin
t0 = celcius_to_kelvin(t0);
t1 = celcius_to_kelvin(t1);

% Calculate constant
Beta_m = (log(r1/r0))/(1/t1 - 1/t0);

% Build anonymous function for resistance 
R = @(T) (r0.*exp(Beta_m .* (1./celcius_to_kelvin(T) - 1/t0)));

% Build anonymous function for temperature
T = @(R) kelvin_to_celcius((1/t0 + log(R./r0)./Beta_m).^(-1));

end

