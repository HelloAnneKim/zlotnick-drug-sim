%% For Calculating K_n_prime
function K_n_prime = Calc_K_n_prime(c)
%Calc_K_n_prime is a function that calculates the K_n_prime value given the
%contacts
%
%   k_4_prime = Calc_K_n_prime(2)
%
%   Refer back to the Zlotnick 1994 paper

delta_G_c_naught = -2.76; % kCal/Mol; the energy of a 2-fold axis
R = 1.987; % 1.987 cal/deg/Mol
T = 298; % K
K_n_prime = exp(-c.*delta_G_c_naught./R./T);
end

