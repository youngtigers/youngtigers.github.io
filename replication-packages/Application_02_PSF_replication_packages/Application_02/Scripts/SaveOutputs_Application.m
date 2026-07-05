% save specified parameters
save Params.mat Params
movefile('Params.mat', 'Results')

% save seives classification
save store_group_sieves.mat store_group_sieves
movefile('store_group_sieves.mat', 'Results')

% save IC 
save store_IC.mat store_IC
movefile('store_IC.mat', 'Results') 

save store_IC1_prime.mat store_IC1_prime
movefile('store_IC1_prime.mat', 'Results')

save store_IC2_prime.mat store_IC2_prime
movefile('store_IC2_prime.mat', 'Results')

save store_Kstar_L1.mat Kstar_L1 store_Kstar_L1
movefile('store_Kstar_L1.mat', 'Results')

% c-u mixture structure under step-5' penalties (1=unique, 2=mixture)
save store_MixD_L1.mat  store_MixD_L1   % c_tilde_lambda = 1 (benchmark)
movefile('store_MixD_L1.mat', 'Results')

save store_MixD_L34.mat store_MixD_L34  % c_tilde_lambda = 3/4
movefile('store_MixD_L34.mat', 'Results')

save store_MixD_L32.mat store_MixD_L32  % c_tilde_lambda = 3/2
movefile('store_MixD_L32.mat', 'Results')

% omegabar: save if it contain values other than NaN 
if ~isempty(store_omegabar_K1(~isnan(store_omegabar_K1)))
    save store_omegabar_K1.mat store_omegabar_K1
    movefile('store_omegabar_K1.mat', 'Results')
elseif ~isempty(store_omegabar_K2(~isnan(store_omegabar_K2)))
    save store_omegabar_K2.mat store_omegabar_K2
    movefile('store_omegabar_K2.mat', 'Results')   
elseif ~isempty(store_omegabar_K3(~isnan(store_omegabar_K3)))
    save store_omegabar_K3.mat store_omegabar_K3
    movefile('store_omegabar_K3.mat', 'Results')
elseif ~isempty(store_omegabar_K4(~isnan(store_omegabar_K4)))
    save store_omegabar_K4.mat store_omegabar_K4
    movefile('store_omegabar_K4.mat', 'Results')
end

% sigvhat: save if it contain values other than NaN
if ~isempty(store_sigvhat_K1(~isnan(store_sigvhat_K1)))
    save store_sigvhat_K1.mat store_sigvhat_K1
    movefile('store_sigvhat_K1.mat', 'Results')
elseif ~isempty(store_sigvhat_K2(~isnan(store_sigvhat_K2)))
    save store_sigvhat_K2.mat store_sigvhat_K2
    movefile('store_sigvhat_K2.mat', 'Results')   
elseif ~isempty(store_sigvhat_K3(~isnan(store_sigvhat_K3)))
    save store_sigvhat_K3.mat store_sigvhat_K3
    movefile('store_sigvhat_K3.mat', 'Results')
elseif ~isempty(store_sigvhat_K4(~isnan(store_sigvhat_K4)))
    save store_sigvhat_K4.mat store_sigvhat_K4
    movefile('store_sigvhat_K4.mat', 'Results')
end
%% frontiers: save the selected-K* store_pi cells (one per group)
if Kstar_L1 == 1
    save store_pi_K1.mat store_pi_K1
    movefile('store_pi_K1.mat', 'Results')
elseif Kstar_L1 == 2
    save store_pi_K2_g1.mat store_pi_K2_g1
    movefile('store_pi_K2_g1.mat', 'Results')
    save store_pi_K2_g2.mat store_pi_K2_g2
    movefile('store_pi_K2_g2.mat', 'Results')
elseif Kstar_L1 == 3
    save store_pi_K3_g1.mat store_pi_K3_g1
    movefile('store_pi_K3_g1.mat', 'Results')
    save store_pi_K3_g2.mat store_pi_K3_g2
    movefile('store_pi_K3_g2.mat', 'Results')
    save store_pi_K3_g3.mat store_pi_K3_g3
    movefile('store_pi_K3_g3.mat', 'Results')
elseif Kstar_L1 == 4
    save store_pi_K4_g1.mat store_pi_K4_g1
    movefile('store_pi_K4_g1.mat', 'Results')
    save store_pi_K4_g2.mat store_pi_K4_g2
    movefile('store_pi_K4_g2.mat', 'Results')
    save store_pi_K4_g3.mat store_pi_K4_g3
    movefile('store_pi_K4_g3.mat', 'Results')
    save store_pi_K4_g4.mat store_pi_K4_g4
    movefile('store_pi_K4_g4.mat', 'Results')
end
