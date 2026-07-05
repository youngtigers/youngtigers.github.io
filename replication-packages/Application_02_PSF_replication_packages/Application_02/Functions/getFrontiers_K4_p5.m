function [hat_alphag1, hat_beta1g1, hat_beta2g1, hat_beta3g1, hat_beta4g1, hat_beta5g1, ...
          hat_alphag2, hat_beta1g2, hat_beta2g2, hat_beta3g2, hat_beta4g2, hat_beta5g2, ...
          hat_alphag3, hat_beta1g3, hat_beta2g3, hat_beta3g3, hat_beta4g3, hat_beta5g3, ...
          hat_alphag4, hat_beta1g4, hat_beta2g4, hat_beta3g4, hat_beta4g4, hat_beta5g4] ...
    = getFrontiers_K4_p5(pi_g1, pi_g2, pi_g3, pi_g4, mLbar_g1, mLbar_g2, mLbar_g3, mLbar_g4, Params)
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

L = mLbar_g2;
B   = sqrt(2) * cos(pi * (t/T) * (1:L-1));
Bmm = [ones(T,1), B];
hat_alphag2 = B   * pi_g2(1:L-1);
hat_beta1g2 = Bmm * pi_g2(L     : 2*L-1);
hat_beta2g2 = Bmm * pi_g2(2*L   : 3*L-1);
hat_beta3g2 = Bmm * pi_g2(3*L   : 4*L-1);
hat_beta4g2 = Bmm * pi_g2(4*L   : 5*L-1);
hat_beta5g2 = Bmm * pi_g2(5*L   : 6*L-1);

L = mLbar_g3;
B   = sqrt(2) * cos(pi * (t/T) * (1:L-1));
Bmm = [ones(T,1), B];
hat_alphag3 = B   * pi_g3(1:L-1);
hat_beta1g3 = Bmm * pi_g3(L     : 2*L-1);
hat_beta2g3 = Bmm * pi_g3(2*L   : 3*L-1);
hat_beta3g3 = Bmm * pi_g3(3*L   : 4*L-1);
hat_beta4g3 = Bmm * pi_g3(4*L   : 5*L-1);
hat_beta5g3 = Bmm * pi_g3(5*L   : 6*L-1);

L = mLbar_g4;
B   = sqrt(2) * cos(pi * (t/T) * (1:L-1));
Bmm = [ones(T,1), B];
hat_alphag4 = B   * pi_g4(1:L-1);
hat_beta1g4 = Bmm * pi_g4(L     : 2*L-1);
hat_beta2g4 = Bmm * pi_g4(2*L   : 3*L-1);
hat_beta3g4 = Bmm * pi_g4(3*L   : 4*L-1);
hat_beta4g4 = Bmm * pi_g4(4*L   : 5*L-1);
hat_beta5g4 = Bmm * pi_g4(5*L   : 6*L-1);
end
