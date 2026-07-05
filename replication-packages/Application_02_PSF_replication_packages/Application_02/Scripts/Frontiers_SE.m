
% storage
hat_alphag1_low = NaN(T,1);
hat_alphag1_high = NaN(T,1);
hat_beta1g1_low = NaN(T,1);
hat_beta1g1_high = NaN(T,1);
hat_beta2g1_low = NaN(T,1);
hat_beta2g1_high = NaN(T,1);
hat_beta3g1_low = NaN(T,1); 
hat_beta3g1_high = NaN(T,1);
hat_beta4g1_low = NaN(T,1);
hat_beta4g1_high = NaN(T,1);
hat_beta5g1_low = NaN(T,1);
hat_beta5g1_high = NaN(T,1);

hat_alphag2_low = NaN(T,1);
hat_alphag2_high = NaN(T,1);
hat_beta1g2_low = NaN(T,1);
hat_beta1g2_high = NaN(T,1);
hat_beta2g2_low = NaN(T,1);
hat_beta2g2_high = NaN(T,1);
hat_beta3g2_low = NaN(T,1); 
hat_beta3g2_high = NaN(T,1);
hat_beta4g2_low = NaN(T,1);
hat_beta4g2_high = NaN(T,1);
hat_beta5g2_low = NaN(T,1);
hat_beta5g2_high = NaN(T,1);

hat_alphag3_low = NaN(T,1);
hat_alphag3_high = NaN(T,1);
hat_beta1g3_low = NaN(T,1);
hat_beta1g3_high = NaN(T,1);
hat_beta2g3_low = NaN(T,1);
hat_beta2g3_high = NaN(T,1);
hat_beta3g3_low = NaN(T,1); 
hat_beta3g3_high = NaN(T,1);
hat_beta4g3_low = NaN(T,1);
hat_beta4g3_high = NaN(T,1);
hat_beta5g3_low = NaN(T,1);
hat_beta5g3_high = NaN(T,1);

hat_alphag4_low = NaN(T,1);
hat_alphag4_high = NaN(T,1);
hat_beta1g4_low = NaN(T,1);
hat_beta1g4_high = NaN(T,1);
hat_beta2g4_low = NaN(T,1);
hat_beta2g4_high = NaN(T,1);
hat_beta3g4_low = NaN(T,1); 
hat_beta3g4_high = NaN(T,1);
hat_beta4g4_low = NaN(T,1);
hat_beta4g4_high = NaN(T,1);
hat_beta5g4_low = NaN(T,1);
hat_beta5g4_high = NaN(T,1);

for t = 1:T
    if Kstar_L1 == 2
    % group 1
    hat_alphag1_low(t,1) = hat_alphag1(t,1) - 1.96*sqrt(S_K2_g1(1,1,t))/sqrt(N_K2_g1*T/mLbar_K2g1);
    hat_alphag1_high(t,1) = hat_alphag1(t,1) + 1.96*sqrt(S_K2_g1(1,1,t))/sqrt(N_K2_g1*T/mLbar_K2g1);
    
    hat_beta1g1_low(t,1) = hat_beta1g1(t,1) - 1.96*sqrt(S_K2_g1(2,2,t))/sqrt(N_K2_g1*T/mLbar_K2g1);
    hat_beta1g1_high(t,1) = hat_beta1g1(t,1) + 1.96*sqrt(S_K2_g1(2,2,t))/sqrt(N_K2_g1*T/mLbar_K2g1);

    hat_beta2g1_low(t,1) = hat_beta2g1(t,1) - 1.96*sqrt(S_K2_g1(3,3,t))/sqrt(N_K2_g1*T/mLbar_K2g1);
    hat_beta2g1_high(t,1) = hat_beta2g1(t,1) + 1.96*sqrt(S_K2_g1(3,3,t))/sqrt(N_K2_g1*T/mLbar_K2g1);
    
    hat_beta3g1_low(t,1) = hat_beta3g1(t,1) - 1.96*sqrt(S_K2_g1(4,4,t))/sqrt(N_K2_g1*T/mLbar_K2g1);
    hat_beta3g1_high(t,1) = hat_beta3g1(t,1) + 1.96*sqrt(S_K2_g1(4,4,t))/sqrt(N_K2_g1*T/mLbar_K2g1);

    hat_beta4g1_low(t,1) = hat_beta4g1(t,1) - 1.96*sqrt(S_K2_g1(5,5,t))/sqrt(N_K2_g1*T/mLbar_K2g1);
    hat_beta4g1_high(t,1) = hat_beta4g1(t,1) + 1.96*sqrt(S_K2_g1(5,5,t))/sqrt(N_K2_g1*T/mLbar_K2g1);
    
    hat_beta5g1_low(t,1) = hat_beta5g1(t,1) - 1.96*sqrt(S_K2_g1(6,6,t))/sqrt(N_K2_g1*T/mLbar_K2g1);
    hat_beta5g1_high(t,1) = hat_beta5g1(t,1) + 1.96*sqrt(S_K2_g1(6,6,t))/sqrt(N_K2_g1*T/mLbar_K2g1);

    % group 2
    hat_alphag2_low(t,1) = hat_alphag2(t,1) - 1.96*sqrt(S_K2_g2(1,1,t))/sqrt(N_K2_g2*T/mLbar_K2g2);
    hat_alphag2_high(t,1) = hat_alphag2(t,1) + 1.96*sqrt(S_K2_g2(1,1,t))/sqrt(N_K2_g2*T/mLbar_K2g2);

    hat_beta1g2_low(t,1) = hat_beta1g2(t,1) - 1.96*sqrt(S_K2_g2(2,2,t))/sqrt(N_K2_g2*T/mLbar_K2g2);
    hat_beta1g2_high(t,1) = hat_beta1g2(t,1) + 1.96*sqrt(S_K2_g2(2,2,t))/sqrt(N_K2_g2*T/mLbar_K2g2);
    
    hat_beta2g2_low(t,1) = hat_beta2g2(t,1) - 1.96*sqrt(S_K2_g2(3,3,t))/sqrt(N_K2_g2*T/mLbar_K2g2);
    hat_beta2g2_high(t,1) = hat_beta2g2(t,1) + 1.96*sqrt(S_K2_g2(3,3,t))/sqrt(N_K2_g2*T/mLbar_K2g2);
    
    hat_beta3g2_low(t,1) = hat_beta3g2(t,1) - 1.96*sqrt(S_K2_g2(4,4,t))/sqrt(N_K2_g2*T/mLbar_K2g2);
    hat_beta3g2_high(t,1) = hat_beta3g2(t,1) + 1.96*sqrt(S_K2_g2(4,4,t))/sqrt(N_K2_g2*T/mLbar_K2g2);

    hat_beta4g2_low(t,1) = hat_beta4g2(t,1) - 1.96*sqrt(S_K2_g2(5,5,t))/sqrt(N_K2_g2*T/mLbar_K2g2);
    hat_beta4g2_high(t,1) = hat_beta4g2(t,1) + 1.96*sqrt(S_K2_g2(5,5,t))/sqrt(N_K2_g2*T/mLbar_K2g2);
    
    hat_beta5g2_low(t,1) = hat_beta5g2(t,1) - 1.96*sqrt(S_K2_g2(6,6,t))/sqrt(N_K2_g2*T/mLbar_K2g2);
    hat_beta5g2_high(t,1) = hat_beta5g2(t,1) + 1.96*sqrt(S_K2_g2(6,6,t))/sqrt(N_K2_g2*T/mLbar_K2g2);

    elseif Kstar_L1 == 3
    % group 1
    hat_alphag1_low(t,1) = hat_alphag1(t,1) - 1.96*sqrt(S_K3_g1(1,1,t))/sqrt(N_K3_g1*T/mLbar_K3g1);
    hat_alphag1_high(t,1) = hat_alphag1(t,1) + 1.96*sqrt(S_K3_g1(1,1,t))/sqrt(N_K3_g1*T/mLbar_K3g1);
    
    hat_beta1g1_low(t,1) = hat_beta1g1(t,1) - 1.96*sqrt(S_K3_g1(2,2,t))/sqrt(N_K3_g1*T/mLbar_K3g1);
    hat_beta1g1_high(t,1) = hat_beta1g1(t,1) + 1.96*sqrt(S_K3_g1(2,2,t))/sqrt(N_K3_g1*T/mLbar_K3g1);

    hat_beta2g1_low(t,1) = hat_beta2g1(t,1) - 1.96*sqrt(S_K3_g1(3,3,t))/sqrt(N_K3_g1*T/mLbar_K3g1);
    hat_beta2g1_high(t,1) = hat_beta2g1(t,1) + 1.96*sqrt(S_K3_g1(3,3,t))/sqrt(N_K3_g1*T/mLbar_K3g1);
    
    hat_beta3g1_low(t,1) = hat_beta3g1(t,1) - 1.96*sqrt(S_K3_g1(4,4,t))/sqrt(N_K3_g1*T/mLbar_K3g1);
    hat_beta3g1_high(t,1) = hat_beta3g1(t,1) + 1.96*sqrt(S_K3_g1(4,4,t))/sqrt(N_K3_g1*T/mLbar_K3g1);

    hat_beta4g1_low(t,1) = hat_beta4g1(t,1) - 1.96*sqrt(S_K3_g1(5,5,t))/sqrt(N_K3_g1*T/mLbar_K3g1);
    hat_beta4g1_high(t,1) = hat_beta4g1(t,1) + 1.96*sqrt(S_K3_g1(5,5,t))/sqrt(N_K3_g1*T/mLbar_K3g1);

    hat_beta5g1_low(t,1) = hat_beta5g1(t,1) - 1.96*sqrt(S_K3_g1(6,6,t))/sqrt(N_K3_g1*T/mLbar_K3g1);
    hat_beta5g1_high(t,1) = hat_beta5g1(t,1) + 1.96*sqrt(S_K3_g1(6,6,t))/sqrt(N_K3_g1*T/mLbar_K3g1);

    % group 2
    hat_alphag2_low(t,1) = hat_alphag2(t,1) - 1.96*sqrt(S_K3_g2(1,1,t))/sqrt(N_K3_g2*T/mLbar_K3g2);
    hat_alphag2_high(t,1) = hat_alphag2(t,1) + 1.96*sqrt(S_K3_g2(1,1,t))/sqrt(N_K3_g2*T/mLbar_K3g2);

    hat_beta1g2_low(t,1) = hat_beta1g2(t,1) - 1.96*sqrt(S_K3_g2(2,2,t))/sqrt(N_K3_g2*T/mLbar_K3g2);
    hat_beta1g2_high(t,1) = hat_beta1g2(t,1) + 1.96*sqrt(S_K3_g2(2,2,t))/sqrt(N_K3_g2*T/mLbar_K3g2);
    
    hat_beta2g2_low(t,1) = hat_beta2g2(t,1) - 1.96*sqrt(S_K3_g2(3,3,t))/sqrt(N_K3_g2*T/mLbar_K3g2);
    hat_beta2g2_high(t,1) = hat_beta2g2(t,1) + 1.96*sqrt(S_K3_g2(3,3,t))/sqrt(N_K3_g2*T/mLbar_K3g2);
    
    hat_beta3g2_low(t,1) = hat_beta3g2(t,1) - 1.96*sqrt(S_K3_g2(4,4,t))/sqrt(N_K3_g2*T/mLbar_K3g2);
    hat_beta3g2_high(t,1) = hat_beta3g2(t,1) + 1.96*sqrt(S_K3_g2(4,4,t))/sqrt(N_K3_g2*T/mLbar_K3g2);

    hat_beta4g2_low(t,1) = hat_beta4g2(t,1) - 1.96*sqrt(S_K3_g2(5,5,t))/sqrt(N_K3_g2*T/mLbar_K3g2);
    hat_beta4g2_high(t,1) = hat_beta4g2(t,1) + 1.96*sqrt(S_K3_g2(5,5,t))/sqrt(N_K3_g2*T/mLbar_K3g2);

    hat_beta5g2_low(t,1) = hat_beta5g2(t,1) - 1.96*sqrt(S_K3_g2(6,6,t))/sqrt(N_K3_g2*T/mLbar_K3g2);
    hat_beta5g2_high(t,1) = hat_beta5g2(t,1) + 1.96*sqrt(S_K3_g2(6,6,t))/sqrt(N_K3_g2*T/mLbar_K3g2);

    % group 3
    hat_alphag3_low(t,1) = hat_alphag3(t,1) - 1.96*sqrt(S_K3_g3(1,1,t))/sqrt(N_K3_g3*T/mLbar_K3g3);
    hat_alphag3_high(t,1) = hat_alphag3(t,1) + 1.96*sqrt(S_K3_g3(1,1,t))/sqrt(N_K3_g3*T/mLbar_K3g3);

    hat_beta1g3_low(t,1) = hat_beta1g3(t,1) - 1.96*sqrt(S_K3_g3(2,2,t))/sqrt(N_K3_g3*T/mLbar_K3g3);
    hat_beta1g3_high(t,1) = hat_beta1g3(t,1) + 1.96*sqrt(S_K3_g3(2,2,t))/sqrt(N_K3_g3*T/mLbar_K3g3);
    
    hat_beta2g3_low(t,1) = hat_beta2g3(t,1) - 1.96*sqrt(S_K3_g3(3,3,t))/sqrt(N_K3_g3*T/mLbar_K3g3);
    hat_beta2g3_high(t,1) = hat_beta2g3(t,1) + 1.96*sqrt(S_K3_g3(3,3,t))/sqrt(N_K3_g3*T/mLbar_K3g3);
    
    hat_beta3g3_low(t,1) = hat_beta3g3(t,1) - 1.96*sqrt(S_K3_g3(4,4,t))/sqrt(N_K3_g3*T/mLbar_K3g3);
    hat_beta3g3_high(t,1) = hat_beta3g3(t,1) + 1.96*sqrt(S_K3_g3(4,4,t))/sqrt(N_K3_g3*T/mLbar_K3g3);

    hat_beta4g3_low(t,1) = hat_beta4g3(t,1) - 1.96*sqrt(S_K3_g3(5,5,t))/sqrt(N_K3_g3*T/mLbar_K3g3);
    hat_beta4g3_high(t,1) = hat_beta4g3(t,1) + 1.96*sqrt(S_K3_g3(5,5,t))/sqrt(N_K3_g3*T/mLbar_K3g3);

    hat_beta5g3_low(t,1) = hat_beta5g3(t,1) - 1.96*sqrt(S_K3_g3(6,6,t))/sqrt(N_K3_g3*T/mLbar_K3g3);
    hat_beta5g3_high(t,1) = hat_beta5g3(t,1) + 1.96*sqrt(S_K3_g3(6,6,t))/sqrt(N_K3_g3*T/mLbar_K3g3);

    elseif Kstar_L1 == 4
    % group 1
    hat_alphag1_low(t,1) = hat_alphag1(t,1) - 1.96*sqrt(S_K4_g1(1,1,t))/sqrt(N_K4_g1*T/mLbar_K4g1);
    hat_alphag1_high(t,1) = hat_alphag1(t,1) + 1.96*sqrt(S_K4_g1(1,1,t))/sqrt(N_K4_g1*T/mLbar_K4g1);
    
    hat_beta1g1_low(t,1) = hat_beta1g1(t,1) - 1.96*sqrt(S_K4_g1(2,2,t))/sqrt(N_K4_g1*T/mLbar_K4g1);
    hat_beta1g1_high(t,1) = hat_beta1g1(t,1) + 1.96*sqrt(S_K4_g1(2,2,t))/sqrt(N_K4_g1*T/mLbar_K4g1);

    hat_beta2g1_low(t,1) = hat_beta2g1(t,1) - 1.96*sqrt(S_K4_g1(3,3,t))/sqrt(N_K4_g1*T/mLbar_K4g1);
    hat_beta2g1_high(t,1) = hat_beta2g1(t,1) + 1.96*sqrt(S_K4_g1(3,3,t))/sqrt(N_K4_g1*T/mLbar_K4g1);
    
    hat_beta3g1_low(t,1) = hat_beta3g1(t,1) - 1.96*sqrt(S_K4_g1(4,4,t))/sqrt(N_K4_g1*T/mLbar_K4g1);
    hat_beta3g1_high(t,1) = hat_beta3g1(t,1) + 1.96*sqrt(S_K4_g1(4,4,t))/sqrt(N_K4_g1*T/mLbar_K4g1);

    hat_beta4g1_low(t,1) = hat_beta4g1(t,1) - 1.96*sqrt(S_K4_g1(5,5,t))/sqrt(N_K4_g1*T/mLbar_K4g1);
    hat_beta4g1_high(t,1) = hat_beta4g1(t,1) + 1.96*sqrt(S_K4_g1(5,5,t))/sqrt(N_K4_g1*T/mLbar_K4g1);
    
    hat_beta5g1_low(t,1) = hat_beta5g1(t,1) - 1.96*sqrt(S_K4_g1(6,6,t))/sqrt(N_K4_g1*T/mLbar_K4g1);
    hat_beta5g1_high(t,1) = hat_beta5g1(t,1) + 1.96*sqrt(S_K4_g1(6,6,t))/sqrt(N_K4_g1*T/mLbar_K4g1);

    % group 2
    hat_alphag2_low(t,1) = hat_alphag2(t,1) - 1.96*sqrt(S_K4_g2(1,1,t))/sqrt(N_K4_g2*T/mLbar_K4g2);
    hat_alphag2_high(t,1) = hat_alphag2(t,1) + 1.96*sqrt(S_K4_g2(1,1,t))/sqrt(N_K4_g2*T/mLbar_K4g2);

    hat_beta1g2_low(t,1) = hat_beta1g2(t,1) - 1.96*sqrt(S_K4_g2(2,2,t))/sqrt(N_K4_g2*T/mLbar_K4g2);
    hat_beta1g2_high(t,1) = hat_beta1g2(t,1) + 1.96*sqrt(S_K4_g2(2,2,t))/sqrt(N_K4_g2*T/mLbar_K4g2);
    
    hat_beta2g2_low(t,1) = hat_beta2g2(t,1) - 1.96*sqrt(S_K4_g2(3,3,t))/sqrt(N_K4_g2*T/mLbar_K4g2);
    hat_beta2g2_high(t,1) = hat_beta2g2(t,1) + 1.96*sqrt(S_K4_g2(3,3,t))/sqrt(N_K4_g2*T/mLbar_K4g2);
    
    hat_beta3g2_low(t,1) = hat_beta3g2(t,1) - 1.96*sqrt(S_K4_g2(4,4,t))/sqrt(N_K4_g2*T/mLbar_K4g2);
    hat_beta3g2_high(t,1) = hat_beta3g2(t,1) + 1.96*sqrt(S_K4_g2(4,4,t))/sqrt(N_K4_g2*T/mLbar_K4g2);

    hat_beta4g2_low(t,1) = hat_beta4g2(t,1) - 1.96*sqrt(S_K4_g2(5,5,t))/sqrt(N_K4_g2*T/mLbar_K4g2);
    hat_beta4g2_high(t,1) = hat_beta4g2(t,1) + 1.96*sqrt(S_K4_g2(5,5,t))/sqrt(N_K4_g2*T/mLbar_K4g2);

    hat_beta5g2_low(t,1) = hat_beta5g2(t,1) - 1.96*sqrt(S_K4_g2(6,6,t))/sqrt(N_K4_g2*T/mLbar_K4g2);
    hat_beta5g2_high(t,1) = hat_beta5g2(t,1) + 1.96*sqrt(S_K4_g2(6,6,t))/sqrt(N_K4_g2*T/mLbar_K4g2);

    % group 3
    hat_alphag3_low(t,1) = hat_alphag3(t,1) - 1.96*sqrt(S_K4_g3(1,1,t))/sqrt(N_K4_g3*T/mLbar_K4g3);
    hat_alphag3_high(t,1) = hat_alphag3(t,1) + 1.96*sqrt(S_K4_g3(1,1,t))/sqrt(N_K4_g3*T/mLbar_K4g3);

    hat_beta1g3_low(t,1) = hat_beta1g3(t,1) - 1.96*sqrt(S_K4_g3(2,2,t))/sqrt(N_K4_g3*T/mLbar_K4g3);
    hat_beta1g3_high(t,1) = hat_beta1g3(t,1) + 1.96*sqrt(S_K4_g3(2,2,t))/sqrt(N_K4_g3*T/mLbar_K4g3);
    
    hat_beta2g3_low(t,1) = hat_beta2g3(t,1) - 1.96*sqrt(S_K4_g3(3,3,t))/sqrt(N_K4_g3*T/mLbar_K4g3);
    hat_beta2g3_high(t,1) = hat_beta2g3(t,1) + 1.96*sqrt(S_K4_g3(3,3,t))/sqrt(N_K4_g3*T/mLbar_K4g3);
    
    hat_beta3g3_low(t,1) = hat_beta3g3(t,1) - 1.96*sqrt(S_K4_g3(4,4,t))/sqrt(N_K4_g3*T/mLbar_K4g3);
    hat_beta3g3_high(t,1) = hat_beta3g3(t,1) + 1.96*sqrt(S_K4_g3(4,4,t))/sqrt(N_K4_g3*T/mLbar_K4g3);

    hat_beta4g3_low(t,1) = hat_beta4g3(t,1) - 1.96*sqrt(S_K4_g3(5,5,t))/sqrt(N_K4_g3*T/mLbar_K4g3);
    hat_beta4g3_high(t,1) = hat_beta4g3(t,1) + 1.96*sqrt(S_K4_g3(5,5,t))/sqrt(N_K4_g3*T/mLbar_K4g3);
    
    hat_beta5g3_low(t,1) = hat_beta5g3(t,1) - 1.96*sqrt(S_K4_g3(6,6,t))/sqrt(N_K4_g3*T/mLbar_K4g3);
    hat_beta5g3_high(t,1) = hat_beta5g3(t,1) + 1.96*sqrt(S_K4_g3(6,6,t))/sqrt(N_K4_g3*T/mLbar_K4g3);
    % group 4
    hat_alphag4_low(t,1) = hat_alphag4(t,1) - 1.96*sqrt(S_K4_g4(1,1,t))/sqrt(N_K4_g4*T/mLbar_K4g4);
    hat_alphag4_high(t,1) = hat_alphag4(t,1) + 1.96*sqrt(S_K4_g4(1,1,t))/sqrt(N_K4_g4*T/mLbar_K4g4);

    hat_beta1g4_low(t,1) = hat_beta1g4(t,1) - 1.96*sqrt(S_K4_g4(2,2,t))/sqrt(N_K4_g4*T/mLbar_K4g4);
    hat_beta1g4_high(t,1) = hat_beta1g4(t,1) + 1.96*sqrt(S_K4_g4(2,2,t))/sqrt(N_K4_g4*T/mLbar_K4g4);
    
    hat_beta2g4_low(t,1) = hat_beta2g4(t,1) - 1.96*sqrt(S_K4_g4(3,3,t))/sqrt(N_K4_g4*T/mLbar_K4g4);
    hat_beta2g4_high(t,1) = hat_beta2g4(t,1) + 1.96*sqrt(S_K4_g4(3,3,t))/sqrt(N_K4_g4*T/mLbar_K4g4);
    
    hat_beta3g4_low(t,1) = hat_beta3g4(t,1) - 1.96*sqrt(S_K4_g4(4,4,t))/sqrt(N_K4_g4*T/mLbar_K4g4);
    hat_beta3g4_high(t,1) = hat_beta3g4(t,1) + 1.96*sqrt(S_K4_g4(4,4,t))/sqrt(N_K4_g4*T/mLbar_K4g4);

    hat_beta4g4_low(t,1) = hat_beta4g4(t,1) - 1.96*sqrt(S_K4_g4(5,5,t))/sqrt(N_K4_g4*T/mLbar_K4g4);
    hat_beta4g4_high(t,1) = hat_beta4g4(t,1) + 1.96*sqrt(S_K4_g4(5,5,t))/sqrt(N_K4_g4*T/mLbar_K4g4);
    
    hat_beta5g4_low(t,1) = hat_beta5g4(t,1) - 1.96*sqrt(S_K4_g4(6,6,t))/sqrt(N_K4_g4*T/mLbar_K4g4);
    hat_beta5g4_high(t,1) = hat_beta5g4(t,1) + 1.96*sqrt(S_K4_g4(6,6,t))/sqrt(N_K4_g4*T/mLbar_K4g4);

    end
end