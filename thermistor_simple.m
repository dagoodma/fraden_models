function [ R T ] = thermistor_simple( r1, t1, r2, t2 )
%THERMISTOR_SIMPLE Simple model for an NTC thermistor.
%   Given two measured resistances, r1 and r2, and known temperatures in
%   celcius t1 and t2, this function will output two anonymous functions
%   that model the resitance as a function of temperature and visa versa.
%
if (r1 == 0 || t1 == 0 || r2 == 0 || t2 == 0)
    error('All arguments must be non-zero values.');
end

% Celcius to kelvin conversion functions
celcius_to_kelvin = @(TC) TC + 273.15;
kelvin_to_celcius = @(TK) TK - 273.15;

% Convert calibration points from celcius to kelvin
t1 = celcius_to_kelvin(t1);
t2 = celcius_to_kelvin(t2);

% Calculate constant
Beta_m = (log(r2/r1))/(1/t2 - 1/t1);

% Build anonymous function for resistance 
R = @(T) (r1.*exp(Beta_m .* (1./celcius_to_kelvin(T) - 1/t1)));

% Build anonymous function for temperature
T = @(R) kelvin_to_celcius((1/t1 + log(R./r1)./Beta_m).^(-1));

end

