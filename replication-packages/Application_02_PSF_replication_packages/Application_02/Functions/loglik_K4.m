%%
%%  L = loglik_K4(omegabar, y, z, sigvhat, pi) 
%%
%%  "loglik_K4" computes negative of log-likelihood function
%%  used in step 4 of the estimation alogorithm for the case where
%%  c - u_i has NO group structure, assuming there's 4 groups for
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


function L = loglik_K4(omegabar, y, z1, z2, z3, z4, ...
    sigvhat_g1, sigvhat_g2, sigvhat_g3, sigvhat_g4, pi_g1, pi_g2, pi_g3, pi_g4, K_sieves ) 

% sort response and covariates into two groups based on HAC
y_g1  =  y(K_sieves==1,:);
[nn1,T1] = size(y_g1); % get the sample size
        
y_g2  =  y(K_sieves==2,:);
[nn2,T2] = size(y_g2); % get the sample size

y_g3  =  y(K_sieves==3,:);
[nn3,T3] = size(y_g3); % get the sample size
        
y_g4  =  y(K_sieves==4,:);
[nn4,T4] = size(y_g4); % get the sample size

% Parameters

epsilon_it1 = NaN(nn1,T1);
epsilon_it2 = NaN(nn2,T2);
epsilon_it3 = NaN(nn3,T3);
epsilon_it4 = NaN(nn4,T4);

% Parameters

c = omegabar(1);
sqrt_s_u = omegabar(2); % sqrt_s_u is the square root of sigma_u; this ensures sigma_u is always positive
s_u = sqrt_s_u^2;   % square it to make it sigma_u

sqrt_s_v1 = sqrt(sigvhat_g1); % sqrt_s_v1 is the square root of sigma_v; note sigvhat = sqrt(sigv2hat_g1)
sqrt_s_v2 = sqrt(sigvhat_g2); % sqrt_s_v2 is the square root of sigma_v; note sigvhat = sqrt(sigv2hat_g2)
sqrt_s_v3 = sqrt(sigvhat_g3); % sqrt_s_v3 is the square root of sigma_v; note sigvhat = sqrt(sigv2hat_g3)
sqrt_s_v4 = sqrt(sigvhat_g4); % sqrt_s_v3 is the square root of sigma_v; note sigvhat = sqrt(sigv2hat_g3)
s_v1 = sqrt_s_v1^2;
s_v2 = sqrt_s_v2^2;
s_v3 = sqrt_s_v3^2;
s_v4 = sqrt_s_v4^2;

% Group 1
s_i1 = sqrt(s_v1^2 + T1*s_u^2);       
% n X T the residuls
for i = 1:nn1
    for t = 1:T1
        epsilon_it1(i,t) = y_g1(i,t) - c - z1(t,:,i)*pi_g1;
    end
end
mu_star_i1 = -s_u^2/s_i1^2*sum(epsilon_it1,2); % n X 1 summation over t
s_star_i1 = sqrt(s_u^2*s_v1^2/s_i1^2);

% Group 2
s_i2 = sqrt(s_v2^2 + T2*s_u^2);       
% n X T the residuls
for i = 1:nn2
    for t = 1:T2
        epsilon_it2(i,t) = y_g2(i,t) - c - z2(t,:,i)*pi_g2;
    end
end
mu_star_i2 = -s_u^2/s_i2^2*sum(epsilon_it2,2); % n X 1 summation over t
s_star_i2 = sqrt(s_u^2*s_v2^2/s_i2^2);

% Group 3
s_i3 = sqrt(s_v3^2 + T3*s_u^2);       
% n X T the residuls
for i = 1:nn3
    for t = 1:T3
        epsilon_it3(i,t) = y_g3(i,t) - c - z3(t,:,i)*pi_g3;
    end
end
mu_star_i3 = -s_u^2/s_i3^2*sum(epsilon_it3,2); % n X 1 summation over t
s_star_i3 = sqrt(s_u^2*s_v3^2/s_i3^2);

% Group 4
s_i4 = sqrt(s_v4^2 + T4*s_u^2);       
% n X T the residuls
for i = 1:nn4
    for t = 1:T4
        epsilon_it4(i,t) = y_g4(i,t) - c - z4(t,:,i)*pi_g4;
    end
end
mu_star_i4 = -s_u^2/s_i4^2*sum(epsilon_it4,2); % n X 1 summation over t
s_star_i4 = sqrt(s_u^2*s_v4^2/s_i4^2);

% part 1 of the loglikelihood
P1_g1 = -nn1*(T1-1)/2*log(s_v1^2) ;
P1_g2 = -nn2*(T2-1)/2*log(s_v2^2) ;
P1_g3 = -nn3*(T3-1)/2*log(s_v3^2) ;
P1_g4 = -nn4*(T4-1)/2*log(s_v4^2) ;

% part 2 
P2_g1 = -(nn1/2)*log(s_i1^2);
P2_g2 = -(nn2/2)*log(s_i2^2);
P2_g3 = -(nn3/2)*log(s_i3^2);
P2_g4 = -(nn4/2)*log(s_i4^2);

% part 3 
P3_g1 = sum(log(1-normcdf(-mu_star_i1/s_star_i1)));
P3_g2 = sum(log(1-normcdf(-mu_star_i2/s_star_i2)));
P3_g3 = sum(log(1-normcdf(-mu_star_i3/s_star_i3)));
P3_g4 = sum(log(1-normcdf(-mu_star_i4/s_star_i4)));

% part 4 
P4_g1 = (1/2)*sum((mu_star_i1/s_star_i1).^2);
P4_g2 = (1/2)*sum((mu_star_i2/s_star_i2).^2);
P4_g3 = (1/2)*sum((mu_star_i3/s_star_i3).^2);
P4_g4 = (1/2)*sum((mu_star_i4/s_star_i4).^2);

% part 5
P5_g1 = - 1/(2*s_v1^2)*sum(sum(epsilon_it1.^2));
P5_g2 = - 1/(2*s_v2^2)*sum(sum(epsilon_it2.^2));
P5_g3 = - 1/(2*s_v3^2)*sum(sum(epsilon_it3.^2));
P5_g4 = - 1/(2*s_v4^2)*sum(sum(epsilon_it4.^2));


% the minus log likelihood function
L1 = -(P1_g1 + P2_g1 + P3_g1 + P4_g1 + P5_g1);
L2 = -(P1_g2 + P2_g2 + P3_g2 + P4_g2 + P5_g2);
L3 = -(P1_g3 + P2_g3 + P3_g3 + P4_g3 + P5_g3);
L4 = -(P1_g4 + P2_g4 + P3_g4 + P4_g4 + P5_g4);

L = L1 + L2 + L3 + L4;
