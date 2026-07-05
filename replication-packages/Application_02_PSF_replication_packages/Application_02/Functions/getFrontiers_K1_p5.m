function [hat_alphag1, hat_beta1g1, hat_beta2g1, hat_beta3g1, hat_beta4g1, hat_beta5g1] ...
    = getFrontiers_K1_p5(pi_g1, mLbar_g1, Params)
T = Params.T;
t = (1:T)';
L = mLbar_g1;
B   = sqrt(2) * cos(pi * (t/T) * (1:L-1));
Bmm = [ones(T,1), B];
hat_alphag1 = B   * pi_g1(1:L-1);
hat_beta1g1 = Bmm * pi_g1(L     : 2*L-1);
hat_beta2g1 = Bmm * pi_g1(2*L   : 3*L-1);
hat_beta3g1 = Bmm * pi_g1(3*L   : 4*L-1);
hat_beta4g1 = Bmm * pi_g1(4*L   : 5*L-1);
hat_beta5g1 = Bmm * pi_g1(5*L   : 6*L-1);
end
