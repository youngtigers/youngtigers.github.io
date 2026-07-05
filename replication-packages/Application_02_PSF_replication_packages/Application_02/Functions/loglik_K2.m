%%
%%  L = loglik_K2(omegabar, y, z, sigvhat, pi) 
%%
%%  "loglik_K2" computes negative of log-likelihood function
%%  used in step 4 of the estimation alogorithm for the case where
%%  c - u_i has NO group structure, assuming there's 2 groups for
%%  vartheta hat
%%
%%  [Model]
%%      y_it = c + alpha_i(t/T) + x_it' beta_i(t/T) + v_it - u_i
%%      
%%      y_it:           scalar
%%      c:              intercept term
%%      alpha_i(t/T):   unknown function
%%      x_it:           p*1 vector
%%      beta_i,l(t/T):  unknown function
%%
%%      v_it ~ N(0, sigv2), iid across t
%%      u_i  ~ |N(0,sigu2)|
%%
%%  [input]
%%      y:      response (N*T matrix)
%%      z:      independent variable ( T*m(p+1)-1*N array ) -1 for intercept
%%      sighat: estimated in step1, scalar
%%      pi:     pi-tilde-hat estimated in step 3 (m(p+1)-1 * 1 vector) -1 for intercept
%%
%%  [output]
%%      L:  negative of logliklihood function evaluated at specfied values


function L = loglik_K2(omegabar, y, z1, z2, sigvhat_g1, sigvhat_g2, pi_g1, pi_g2, K_sieves ) 

% sort response and covariates into two groups based on HAC
y_g1  =  y(K_sieves==1,:);
[nn1,T1] = size(y_g1); % get the sample size
        
y_g2  =  y(K_sieves==2,:);
[nn2,T2] = size(y_g2); % get the sample size

epsilon_it1 = NaN(nn1,T1);
epsilon_it2 = NaN(nn2,T2);

% Parameters
c = omegabar(1);
sqrt_s_u = omegabar(2); % sqrt_s_u is the square root of sigma_u; this ensures sigma_u is always positive
s_u = sqrt_s_u^2;   % square it to make it sigma_u

sqrt_s_v1 = sqrt(sigvhat_g1); % sqrt_s_v1 is the square root of sigma_v; note sigvhat = sqrt(sigv2hat_g1)
sqrt_s_v2 = sqrt(sigvhat_g2); % sqrt_s_v1 is the square root of sigma_v; note sigvhat = sqrt(sigv2hat_g2)
s_v1 = sqrt_s_v1^2;
s_v2 = sqrt_s_v2^2;


% Group 1
s_i1 = sqrt(s_v1^2 + T1*s_u^2);
for i = 1:nn1
    for t = 1:T1
        epsilon_it1(i,t) = y_g1(i,t) - c - z1(t,:,i)*pi_g1;
    end
end
mu_star_i1 = -s_u^2/s_i1^2*sum(epsilon_it1,2); % n X 1 summation over t
s_star_i1 = sqrt(s_u^2*s_v1^2/s_i1^2);

% Group 2
s_i2 = sqrt(s_v2^2 + T2*s_u^2);
for i = 1:nn2
    for t = 1:T2
        epsilon_it2(i,t) = y_g2(i,t) - c - z2(t,:,i)*pi_g2;
    end
end
mu_star_i2 = -s_u^2/s_i2^2*sum(epsilon_it2,2); % n X 1 summation over t
s_star_i2 = sqrt(s_u^2*s_v2^2/s_i2^2);

% part 1 of the loglikelihood
%P1 = -nn*(T-1)/2*log(s_v^2) ;

P1_g1 = -nn1*(T1-1)/2*log(s_v1^2) ;
P1_g2 = -nn2*(T2-1)/2*log(s_v2^2) ;

% part 2 
%P2 = -(nn/2)*log(s_i^2);

P2_g1 = -(nn1/2)*log(s_i1^2);
P2_g2 = -(nn2/2)*log(s_i2^2);

% part 3 
%P3 = sum(log(1-normcdf(-mu_star_i/s_star_i)));

P3_g1 = sum(log(1-normcdf(-mu_star_i1/s_star_i1)));
P3_g2 = sum(log(1-normcdf(-mu_star_i2/s_star_i2)));

% part 4 
P4_g1 = (1/2)*sum((mu_star_i1/s_star_i1).^2);
P4_g2 = (1/2)*sum((mu_star_i2/s_star_i2).^2);

% part 5
%P5 = - 1/(2*s_v^2)*sum(sum(epsilon_it.^2));

P5_g1 = - 1/(2*s_v1^2)*sum(sum(epsilon_it1.^2));
P5_g2 = - 1/(2*s_v2^2)*sum(sum(epsilon_it2.^2));

% the minus log likelihood function
% L = -(P1 + P2 + P3 + P4 + P5);

L1 = -(P1_g1 + P2_g1 + P3_g1 + P4_g1 + P5_g1);
L2 = -(P1_g2 + P2_g2 + P3_g2 + P4_g2 + P5_g2);

L = L1 + L2;
end
