function [ R T ] = thermistor_fraden( r0, t0, r1, t1, ta, ra, tb, rb, tc, rc )
%THERMISTOR_FRADEN Fraden model for an NTC thermistor.
%   Given two measured resistances r0 and r1 at temperatures in
%   celcius t0 and t1, as well three characterizing restances ra, rb, and rc
%   at temperatures ta, tb, and tc in celcius (typically from the devices
%   datasheet), this function will output two anonymous functions
%   that model the resitance as a function of temperature and visa versa.
%
if (r0 == 0 || t0 == 0 || r1 == 0 || t1 == 0)
    error('All arguments must be non-zero values.');
end

% Celcius to kelvin conversion functions
celcius_to_kelvin = @(TC) TC + 273.15;
kelvin_to_celcius = @(TK) TK - 273.15;

% Celcius to kelvin
t0 = celcius_to_kelvin(t0);
t1 = celcius_to_kelvin(t1);
ta = celcius_to_kelvin(ta);
tb = celcius_to_kelvin(tb);
tc = celcius_to_kelvin(tc);

% Calculate constants
Beta_m = (log(r1/r0))/(1/t1 - 1/t0);
Beta_x = (log(rc/rb)/(1/tc - 1/tb));
Beta_y = (log(ra/rb)/(1/ta - 1/tb));
gamma = (Beta_x/Beta_y - 1) * (1/(tc-ta));


% Build anonymous function for resistance 
R = @(T) (r0.*exp(Beta_m .* (1 + gamma .*(celcius_to_kelvin(T) - t0) .* (1./celcius_to_kelvin(T) - 1/t0))));

% Build anonymous function for temperature
% T_r is an inner function for T
T_r = @(R) (1/t0 + log(R./r0)./Beta_m).^(-1);
T=@(R) kelvin_to_celcius((1/t0 + (log(R./r0))./( Beta_m * (1 - gamma.*(t1 - T_r(R))) ) ).^(-1));

end

