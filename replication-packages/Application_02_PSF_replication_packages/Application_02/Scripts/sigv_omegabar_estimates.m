%% sigv and omegabar estimates
if Kstar_L1 == 1 && store_MixD_L1 == 1

% sigv estimate
fprintf('\n sigv estimate \n');
disp([store_sigvhat_K1])

% Convert sqrt(sigu) into sigu
store_omegabar_K1(:,2) = store_omegabar_K1(:,2).^2; % sigu1

% Mean omegabar = [c1, sigu1]'
fprintf('\n omegabar estimate \n');
disp([store_omegabar_K1])

%%
elseif Kstar_L1 == 1 && store_MixD_L1 == 2
% sigv estimate
fprintf('\n sigv estimate \n');
disp([store_sigvhat_K1])

% omegabar_prime estimates
% Convert sqrt(sigu) into sigu
store_omegabar_K1_prime(:,2) = store_omegabar_K1_prime(:,2).^2; % sigu1
store_omegabar_K1_prime(:,4) = store_omegabar_K1_prime(:,4).^2; % sigu2

% omegabar = [c1, sigu1, c2, sigu2, tau]'
fprintf('\n omegabar estimates \n');
disp([store_omegabar_K1_prime])

%%
elseif Kstar_L1 == 2 && store_MixD_L1 == 1

% sigv estimate
fprintf('\n sigv estimate \n');
disp([store_sigvhat_K2])

% Convert sqrt(sigu) into sigu
store_omegabar_K2(:,2) = store_omegabar_K2(:,2).^2; % sigu1

% Mean omegabar = [c1, sigu1]'
fprintf('\n omegabar estimate \n');
disp([store_omegabar_K2])

%%
elseif Kstar_L1 == 2 && store_MixD_L1 == 2

% sigv estimate
fprintf('\n sigv estimate \n');
disp([store_sigvhat_K2])
fprintf('\n se(sigv) estimate \n');
disp([SE_sigvhat_K2_g1, SE_sigvhat_K2_g2])

% Display order: [tau, c1, sigu1, c2, sigu2]
% The optimizer labels group 1/2 in the opposite order to the paper, with tau = 1-paper_tau,
% so we remap: tau=1-pos5, c1=pos3, sigu1=pos4^2, c2=pos1, sigu2=pos2^2
omegabar_K2_prime_disp = [1 - store_omegabar_K2_prime(5), ...
                               store_omegabar_K2_prime(3), ...
                               store_omegabar_K2_prime(4).^2, ...
                               store_omegabar_K2_prime(1), ...
                               store_omegabar_K2_prime(2).^2];

fprintf('\n omegabar estimates: [tau, c1, sigu1, c2, sigu2] \n');
disp(omegabar_K2_prime_disp)

% se(omegabar): [se_tau, se_c1, se_sqrt_sigu1, se_c2, se_sqrt_sigu2]
fprintf('\n se(omegabar estimates) \n');
se_omegabar_prime = se_omegabar_K2_prime(store_omegabar_K2_prime, y, ZimLbar_K2_g1, ZimLbar_K2_g2, ...
        sqrt(sigv2hat_K2_g1), sqrt(sigv2hat_K2_g2), pi_K2_g1, pi_K2_g2, Ksieves2);
disp(se_omegabar_prime)

%%
elseif Kstar_L1 == 3 && store_MixD_L1 == 1

% sigv estimate
fprintf('\n sigv estimate \n');
disp([store_sigvhat_K3])

% Convert sqrt(sigu) into sigu
store_omegabar_K3(:,2) = store_omegabar_K3(:,2).^2; % sigu1

% Mean omegabar = [c1, sigu1]'
fprintf('\n omegabar estimate \n');
disp([store_omegabar_K3])

%%
elseif Kstar_L1 == 3 && store_MixD_L1 == 2

% sigv estimate
fprintf('\n sigv estimate \n');
disp([store_sigvhat_K3])
disp([store_sigvhat_K3(1)./sqrt(sum(Ksieves3(:)==1)*T), ...
      store_sigvhat_K3(2)./sqrt(sum(Ksieves3(:)==2)*T), ...
      store_sigvhat_K3(3)./sqrt(sum(Ksieves3(:)==3)*T)])
     
% omegabar_prime estimates
% Convert sqrt(sigu) into sigu
store_omegabar_K3_prime(:,2) = store_omegabar_K3_prime(:,2).^2; % sigu1
store_omegabar_K3_prime(:,4) = store_omegabar_K3_prime(:,4).^2; % sigu2

% omegabar = [c1, sigu1, c2, sigu2, tau]'
fprintf('\n omegabar estimates \n');
disp([store_omegabar_K3_prime])

%%
elseif Kstar_L1 == 4 && store_MixD_L1 == 1

% sigv estimate
fprintf('\n sigv estimate \n');
disp([store_sigvhat_K4])

% Convert sqrt(sigu) into sigu
store_omegabar_K4(:,2) = store_omegabar_K4(:,2).^2; % sigu1

% Mean omegabar = [c1, sigu1]'
fprintf('\n omegabar estimate \n');
disp([store_omegabar_K4])

%%
elseif Kstar_L1 == 4 && store_MixD_L1 == 2

% sigv estimate
fprintf('\n sigv estimate \n');
disp([store_sigvhat_K4])
    
% omegabar_prime estimates
% Convert sqrt(sigu) into sigu
store_omegabar_K4_prime(:,2) = store_omegabar_K4_prime(:,2).^2; % sigu1
store_omegabar_K4_prime(:,4) = store_omegabar_K4_prime(:,4).^2; % sigu2

% omegabar = [c1, sigu1, c2, sigu2, tau]'
fprintf('\n omegabar estimates \n');
disp([store_omegabar_K4_prime])
end