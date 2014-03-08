function [ R T ] = thermistor_fraden( r1, t1, r2, t2, ra, ta, rb, tb, rc, tc )
%THERMISTOR_FRADEN Fraden model for an NTC thermistor.
%   Given two measured resistances r1 and r2 at temperatures in
%   celcius t1 and t2, as well three characterizing restances ra, rb, and rc
%   at temperatures ta, tb, and tc in celcius (typically from the device's
%   datasheet), this function will output two anonymous functions
%   that model the resitance as a function of temperature and visa versa.
%
if (r1 == 0 || t1 == 0 || r2 == 0 || t2 == 0 ...
    || ra == 0 || ta == 0 || rb == 0 || tb == 0 || rc == 0 || tc == 0)
    error('All arguments must be non-zero values.');
end

% Celcius to kelvin conversion functions
celcius_to_kelvin = @(TC) TC + 273.15;
kelvin_to_celcius = @(TK) TK - 273.15;

% Celcius to kelvin
t1 = celcius_to_kelvin(t1);
t2 = celcius_to_kelvin(t2);
ta = celcius_to_kelvin(ta);
tb = celcius_to_kelvin(tb);
tc = celcius_to_kelvin(tc);

% Calculate constants
Beta_m = (log(r2/r1))/(1/t2 - 1/t1);
Beta_x = (log(rc/rb)/(1/tc - 1/tb));
Beta_y = (log(ra/rb)/(1/ta - 1/tb));
gamma = (Beta_x/Beta_y - 1) * (1/(tc-ta));


% Build anonymous function for resistance 
R = @(T) (r1.*exp(Beta_m .* (1 + gamma .*(celcius_to_kelvin(T) - t1) .* (1./celcius_to_kelvin(T) - 1/t1))));

% Build anonymous function for temperature
% T_r is an inner function for T
T_r = @(R) (1/t1 + log(R./r1)./Beta_m).^(-1);
T=@(R) kelvin_to_celcius((1/t1 + (log(R./r1))./( Beta_m * (1 - gamma.*(t2 - T_r(R))) ) ).^(-1));

end

