%%
filename = 'Application';

if Kstar_L1 == 1
    pi_K1 = store_pi_K1{mLbar_K1};
    [hat_alphag1, hat_beta1g1, hat_beta2g1, hat_beta3g1, hat_beta4g1, hat_beta5g1] ...
        = getFrontiers_K1_p5(pi_K1, mLbar_K1, Params);

    save(['Results/hat_alphag1_' filename '.mat'],'hat_alphag1')
    save(['Results/hat_beta1g1_' filename '.mat'],'hat_beta1g1')
    save(['Results/hat_beta2g1_' filename '.mat'],'hat_beta2g1')
    save(['Results/hat_beta3g1_' filename '.mat'],'hat_beta3g1')
    save(['Results/hat_beta4g1_' filename '.mat'],'hat_beta4g1')
    save(['Results/hat_beta5g1_' filename '.mat'],'hat_beta5g1')

elseif Kstar_L1 == 2
    pi_K2_g1 = store_pi_K2_g1{mLbar_K2g1};
    pi_K2_g2 = store_pi_K2_g2{mLbar_K2g2};
    [hat_alphag1, hat_beta1g1, hat_beta2g1, hat_beta3g1, hat_beta4g1, hat_beta5g1,...
     hat_alphag2, hat_beta1g2, hat_beta2g2, hat_beta3g2, hat_beta4g2, hat_beta5g2] ...
        = getFrontiers_K2_p5(pi_K2_g1, pi_K2_g2, mLbar_K2g1, mLbar_K2g2, Params);

    save(['Results/hat_alphag1_' filename '.mat'],'hat_alphag1')
    save(['Results/hat_beta1g1_' filename '.mat'],'hat_beta1g1')
    save(['Results/hat_beta2g1_' filename '.mat'],'hat_beta2g1')
    save(['Results/hat_beta3g1_' filename '.mat'],'hat_beta3g1')
    save(['Results/hat_beta4g1_' filename '.mat'],'hat_beta4g1')
    save(['Results/hat_beta5g1_' filename '.mat'],'hat_beta5g1')

    save(['Results/hat_alphag2_' filename '.mat'],'hat_alphag2')
    save(['Results/hat_beta1g2_' filename '.mat'],'hat_beta1g2')
    save(['Results/hat_beta2g2_' filename '.mat'],'hat_beta2g2')
    save(['Results/hat_beta3g2_' filename '.mat'],'hat_beta3g2')
    save(['Results/hat_beta4g2_' filename '.mat'],'hat_beta4g2')
    save(['Results/hat_beta5g2_' filename '.mat'],'hat_beta5g2')

elseif Kstar_L1 == 3
    pi_K3_g1 = store_pi_K3_g1{mLbar_K3g1};
    pi_K3_g2 = store_pi_K3_g2{mLbar_K3g2};
    pi_K3_g3 = store_pi_K3_g3{mLbar_K3g3};
    [hat_alphag1, hat_beta1g1, hat_beta2g1, hat_beta3g1, hat_beta4g1, hat_beta5g1, ...
     hat_alphag2, hat_beta1g2, hat_beta2g2, hat_beta3g2, hat_beta4g2, hat_beta5g2, ...
     hat_alphag3, hat_beta1g3, hat_beta2g3, hat_beta3g3, hat_beta4g3, hat_beta5g3] ...
        = getFrontiers_K3_p5(pi_K3_g1, pi_K3_g2, pi_K3_g3, ...
                             mLbar_K3g1, mLbar_K3g2, mLbar_K3g3, Params);

    save(['Results/hat_alphag1_' filename '.mat'],'hat_alphag1')
    save(['Results/hat_beta1g1_' filename '.mat'],'hat_beta1g1')
    save(['Results/hat_beta2g1_' filename '.mat'],'hat_beta2g1')
    save(['Results/hat_beta3g1_' filename '.mat'],'hat_beta3g1')
    save(['Results/hat_beta4g1_' filename '.mat'],'hat_beta4g1')
    save(['Results/hat_beta5g1_' filename '.mat'],'hat_beta5g1')

    save(['Results/hat_alphag2_' filename '.mat'],'hat_alphag2')
    save(['Results/hat_beta1g2_' filename '.mat'],'hat_beta1g2')
    save(['Results/hat_beta2g2_' filename '.mat'],'hat_beta2g2')
    save(['Results/hat_beta3g2_' filename '.mat'],'hat_beta3g2')
    save(['Results/hat_beta4g2_' filename '.mat'],'hat_beta4g2')
    save(['Results/hat_beta5g2_' filename '.mat'],'hat_beta5g2')

    save(['Results/hat_alphag3_' filename '.mat'],'hat_alphag3')
    save(['Results/hat_beta1g3_' filename '.mat'],'hat_beta1g3')
    save(['Results/hat_beta2g3_' filename '.mat'],'hat_beta2g3')
    save(['Results/hat_beta3g3_' filename '.mat'],'hat_beta3g3')
    save(['Results/hat_beta4g3_' filename '.mat'],'hat_beta4g3')
    save(['Results/hat_beta5g3_' filename '.mat'],'hat_beta5g3')

elseif Kstar_L1 == 4
    pi_K4_g1 = store_pi_K4_g1{mLbar_K4g1};
    pi_K4_g2 = store_pi_K4_g2{mLbar_K4g2};
    pi_K4_g3 = store_pi_K4_g3{mLbar_K4g3};
    pi_K4_g4 = store_pi_K4_g4{mLbar_K4g4};
    [hat_alphag1, hat_beta1g1, hat_beta2g1, hat_beta3g1, hat_beta4g1, hat_beta5g1,...
     hat_alphag2, hat_beta1g2, hat_beta2g2, hat_beta3g2, hat_beta4g2, hat_beta5g2,...
     hat_alphag3, hat_beta1g3, hat_beta2g3, hat_beta3g3, hat_beta4g3, hat_beta5g3,...
     hat_alphag4, hat_beta1g4, hat_beta2g4, hat_beta3g4, hat_beta4g4, hat_beta5g4] ...
        = getFrontiers_K4_p5(pi_K4_g1, pi_K4_g2, pi_K4_g3, pi_K4_g4, ...
                             mLbar_K4g1, mLbar_K4g2, mLbar_K4g3, mLbar_K4g4, Params);

    save(['Results/hat_alphag1_' filename '.mat'],'hat_alphag1')
    save(['Results/hat_beta1g1_' filename '.mat'],'hat_beta1g1')
    save(['Results/hat_beta2g1_' filename '.mat'],'hat_beta2g1')
    save(['Results/hat_beta3g1_' filename '.mat'],'hat_beta3g1')
    save(['Results/hat_beta4g1_' filename '.mat'],'hat_beta4g1')
    save(['Results/hat_beta5g1_' filename '.mat'],'hat_beta5g1')

    save(['Results/hat_alphag2_' filename '.mat'],'hat_alphag2')
    save(['Results/hat_beta1g2_' filename '.mat'],'hat_beta1g2')
    save(['Results/hat_beta2g2_' filename '.mat'],'hat_beta2g2')
    save(['Results/hat_beta3g2_' filename '.mat'],'hat_beta3g2')
    save(['Results/hat_beta4g2_' filename '.mat'],'hat_beta4g2')
    save(['Results/hat_beta5g2_' filename '.mat'],'hat_beta5g2')

    save(['Results/hat_alphag3_' filename '.mat'],'hat_alphag3')
    save(['Results/hat_beta1g3_' filename '.mat'],'hat_beta1g3')
    save(['Results/hat_beta2g3_' filename '.mat'],'hat_beta2g3')
    save(['Results/hat_beta3g3_' filename '.mat'],'hat_beta3g3')
    save(['Results/hat_beta4g3_' filename '.mat'],'hat_beta4g3')
    save(['Results/hat_beta5g3_' filename '.mat'],'hat_beta5g3')

    save(['Results/hat_alphag4_' filename '.mat'],'hat_alphag4')
    save(['Results/hat_beta1g4_' filename '.mat'],'hat_beta1g4')
    save(['Results/hat_beta2g4_' filename '.mat'],'hat_beta2g4')
    save(['Results/hat_beta3g4_' filename '.mat'],'hat_beta3g4')
    save(['Results/hat_beta4g4_' filename '.mat'],'hat_beta4g4')
    save(['Results/hat_beta5g4_' filename '.mat'],'hat_beta5g4')
end

% Save outputs
run SaveOutputs_Application
