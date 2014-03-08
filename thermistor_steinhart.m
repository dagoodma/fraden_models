function [ R T ] = thermistor_steinhart( r1, t1, r2, t2, r3, t3 )
%THERMISTOR_STEINHART Steinhart-Hart model for an NTC thermistor.
%   Given three measured resistances, r1, r2, and r3 and known temperatures in
%   celcius t1, t2, and t3, this function will output two anonymous functions
%   that model the resitance as a function of temperature and visa versa.
%
if (r1 == 0 || t1 == 0 || r2 == 0 || t2 == 0 || r3 == 0 || t3 == 0)
    error('All arguments must be non-zero values.');
end

% Celcius to kelvin conversion functions
celcius_to_kelvin = @(TC) TC + 273.15;
kelvin_to_celcius = @(TK) TK - 273.15;

% Convert calibration points from celcius to kelvin
t1 = celcius_to_kelvin(t1);
t2 = celcius_to_kelvin(t2);
t3 = celcius_to_kelvin(t3);

% Constants
Y1=1/t1; Y2=1/t2; Y3=1/t3;
L1 = log(r1); L2 = log(r2); L3 = log(r3);
gamma2  = (Y2 - Y1)/(L2 - L1);
gamma3  = (Y3 - Y1)/(L3 - L1);

C = (gamma3 - gamma2)/(L3 - L2)*(L1 + L2 + L3)^(-1)'
B = gamma2 - C*(L1^2 + L1*L2 + L2^2);
A = Y1 - (B + L1^2*C)*L1;

% Build anonymous function for resistance
y = @(T) (A - (1./T))./(2*C);
x = @(T) sqrt((B/(3*C))^3 + y(T).^2)
R = @(T) exp((x(T) - y(T)).^(1/3) - (x(T) + y(T)).^(1/3));

% Build anonymous function for temperature
T = @(R) kelvin_to_celcius( (A + B.*log(R) + C.*(log(R)).^3).^(-1) );


end

