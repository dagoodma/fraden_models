function [ R T ] = thermistor_steinhart( r0, t0, r1, t1, r2, t2 )
%THERMISTOR_STEINHART Steinhart-Hart model for an NTC thermistor.
%   Given three measured resistances, r0, r1, and r2 and known temperatures in
%   celcius t0, t1, and t2, this function will output two anonymous functions
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
t2 = celcius_to_kelvin(t2);

% Constants
G = t0^(-1) - t1^(-1);
H = t0^(-1) - t2^(-1);
F = log(r0) - log(r2);
Z = log(r0) - log(r1);

C = (G - (Z*H)/F)*((log(r0)^3-log(r1)^3) - (Z/F)*(log(r0)^3 - log(r2)^3))^(-1);
B = (Z^(-1))*(G-C*(log(r0^3)-log(r1^3)));
A = t0^(-1) - C*log(r0)^3-B*log(r0);


% Build anonymous function for resistance
R = @(T) exp(((B.^3./(27.*C.^3) + (A.*celcius_to_kelvin(T) - 1).^2/(4.*C.^2.*celcius_to_kelvin(T).^2)).^(1/2) - (A.*celcius_to_kelvin(T) - 1)./(2*C.*celcius_to_kelvin(T))).^(1/3) - B./(3*C.*((B.^3./(27*C.^3) + (A.*celcius_to_kelvin(T) - 1).^2./(4*C.^2.*celcius_to_kelvin(T).^2)).^(1/2) - (A.*celcius_to_kelvin(T) - 1)./(2*C.*celcius_to_kelvin(T))).^(1/3)));

% Build anonymous function for temperature
T = @(R) kelvin_to_celcius((A + B.*log(R) + C.*(log(R)).^3).^(-1));


end

