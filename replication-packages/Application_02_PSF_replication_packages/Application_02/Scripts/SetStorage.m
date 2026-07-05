%% Set storage

% set number of regressors
p = Params.p;

K_group = Params.K_group; 

% store pihat per (K, group) as cell arrays indexed by m lower bar (mLbar).
% Only the chosen-m slot is filled (in main_application); using cells avoids
% the ~100 NaN flat store_pi_K*_g*_<m> variables the old scheme created.
M = 12;   % max m lower bar
store_pi_K1    = cell(1,M);
store_pi_K2_g1 = cell(1,M);
store_pi_K2_g2 = cell(1,M);
store_pi_K3_g1 = cell(1,M);
store_pi_K3_g2 = cell(1,M);
store_pi_K3_g3 = cell(1,M);
store_pi_K4_g1 = cell(1,M);
store_pi_K4_g2 = cell(1,M);
store_pi_K4_g3 = cell(1,M);
store_pi_K4_g4 = cell(1,M);

% store sigvhat := \hat{\sigma}_{v(k|K)} 
store_sigvhat_K1 = NaN(1,1);
store_sigvhat_K2 = NaN(1,2);
store_sigvhat_K3 = NaN(1,3);
store_sigvhat_K4 = NaN(1,4);

% store omegabar := \hat{\varrho}_{k|K} = [c,sigu]
store_omegabar_K1 = NaN(1,2);
store_omegabar_K2 = NaN(1,2);
store_omegabar_K3 = NaN(1,2);
store_omegabar_K4 = NaN(1,2);

% store omegabar := \hat{\varrho}_{k|K} = [c1,sigu1,c2,sigu2,tau]
store_omegabar_K1_prime = NaN(1,5);
store_omegabar_K2_prime = NaN(1,5);
store_omegabar_K3_prime = NaN(1,5);
store_omegabar_K4_prime = NaN(1,5);

store_vthetai = NaN(N,3);
store_group_sieves = NaN(N,1,K_group);

% store IC
store_IC = NaN(K_group,1);
store_IC1_prime = NaN(K_group,1);
store_IC2_prime = NaN(K_group,1);

store_Kstar = NaN(1,1);
store_MixD  = NaN(1,1);





