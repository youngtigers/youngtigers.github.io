%%
%%  Lprime = loglik_K3_prime(omegabar, y, z, sigvhat, pi) 
%%
%%  "loglik_K1_prime" computes negative of log-likelihood function
%%  used in step 4' of the estimation alogorithm for the case where
%%  c - u_i has a group structure, assuming there's 3 groups for
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
%%      Lprime:  negative of logliklihood function evaluated at specfied values

function Lprime = loglik_K3_prime(omegabar, y, z1, z2, z3, sigvhat_g1, sigvhat_g2, sigvhat_g3, pi_g1, pi_g2, pi_g3, K_sieves) 

% sort response and covariates into two groups based on HAC
y_g1  =  y(K_sieves==1,:);
[nn1,T1] = size(y_g1); % get the sample size
        
y_g2  =  y(K_sieves==2,:);
[nn2,T2] = size(y_g2); % get the sample size

y_g3  =  y(K_sieves==3,:);
[nn3,T3] = size(y_g3); % get the sample size
        
% initialize
epsilon_it11 = NaN(nn1,T1);
epsilon_it12 = NaN(nn1,T1);
epsilon_it21 = NaN(nn2,T2);
epsilon_it22 = NaN(nn2,T2);
epsilon_it31 = NaN(nn3,T3);
epsilon_it32 = NaN(nn3,T3);

% Parameters
c1        = omegabar(1);
sqrt_s_u1 = omegabar(2); % sqrt_s_u is the square root of sigma_u1; this ensures sigma_u1 is always positive
c2        = omegabar(3);
sqrt_s_u2 = omegabar(4); % sqrt_s_u is the square root of sigma_u2; this ensures sigma_u2 is always positive
tau       = omegabar(5);

% Parameters
sqrt_s_v1 = sqrt(sigvhat_g1); % sqrt_s_v1 is the square root of sigma_v; note sigvhat = sqrt(sigv2hat_g1)
sqrt_s_v2 = sqrt(sigvhat_g2); % sqrt_s_v2 is the square root of sigma_v; note sigvhat = sqrt(sigv2hat_g2)
sqrt_s_v3 = sqrt(sigvhat_g3); % sqrt_s_v3 is the square root of sigma_v; note sigvhat = sqrt(sigv2hat_g3)
s_v1 = sqrt_s_v1^2;
s_v2 = sqrt_s_v2^2;
s_v3 = sqrt_s_v3^2;

s_u1 = sqrt_s_u1^2;   % square it to make it sigma_u1
s_u2 = sqrt_s_u2^2;   % square it to make it sigma_u2


% Group 1, a1, u1
s_i11 = sqrt(s_v1^2 + T1*s_u1^2);
for i = 1:nn1
    for t = 1:T1
        epsilon_it11(i,t) = y_g1(i,t) - (c1) - z1(t,:,i)*pi_g1;  % n X T the residuls
    end
end
mu_star_i11 = -s_u1^2/s_i11^2*sum(epsilon_it11,2); % n X 1 summation over t
s_star_i11 = sqrt(s_u1^2*s_v1^2/s_i11^2);

% Group 1, a2, u2
s_i12 = sqrt(s_v1^2 + T1*s_u2^2);
%epsilon_it12 = y_g1 - c2 - x1_g1*b1_g1 - x2_g1*b2_g1;         % n X T the residuls
for i = 1:nn1
    for t = 1:T1
        epsilon_it12(i,t) = y_g1(i,t) - (c2) - z1(t,:,i)*pi_g1;
    end
end
mu_star_i12 = -s_u2^2/s_i12^2*sum(epsilon_it12,2); % n X 1 summation over t
s_star_i12 = sqrt(s_u2^2*s_v1^2/s_i12^2);

% Group 2,  a1, u1
s_i21 = sqrt(s_v2^2 + T2*s_u1^2);
%epsilon_it21 = y_g2 - a1 - x1_g2*b1_g2 - x2_g2*b2_g2;         % n X T the residuls
for i = 1:nn2
    for t = 1:T2
        epsilon_it21(i,t) = y_g2(i,t) - (c1) - z2(t,:,i)*pi_g2;
    end
end
mu_star_i21 = -s_u1^2/s_i21^2*sum(epsilon_it21,2); % n X 1 summation over t
s_star_i21 = sqrt(s_u1^2*s_v2^2/s_i21^2);

% Group 2,  a2, u2
s_i22 = sqrt(s_v2^2 + T2*s_u2^2);
%epsilon_it22 = y_g2 - a2 - x1_g2*b1_g2 - x2_g2*b2_g2;         % n X T the residuls
for i = 1:nn2
    for t = 1:T2
        epsilon_it22(i,t) = y_g2(i,t) - (c2) - z2(t,:,i)*pi_g2;
    end
end
mu_star_i22 = -s_u2^2/s_i22^2*sum(epsilon_it22,2); % n X 1 summation over t
s_star_i22 = sqrt(s_u2^2*s_v2^2/s_i22^2);


% Group 3,  a1, u1
s_i31 = sqrt(s_v3^2 + T3*s_u1^2);
%epsilon_it31 = y_g3 - a1 - x1_g3*b1_g3 - x2_g3*b2_g3;         % n X T the residuls
for i = 1:nn3
    for t = 1:T3
        epsilon_it31(i,t) = y_g3(i,t) - (c1) - z3(t,:,i)*pi_g3;
    end
end
mu_star_i31 = -s_u1^2/s_i31^2*sum(epsilon_it31,2); % n X 1 summation over t
s_star_i31 = sqrt(s_u1^2*s_v3^2/s_i31^2);

% Group 3,  a2, u2
s_i32 = sqrt(s_v3^2 + T3*s_u2^2);
%epsilon_it32 = y_g3 - a2 - x1_g3*b1_g3 - x2_g3*b2_g3;         % n X T the residuls
for i = 1:nn3
    for t = 1:T3
        epsilon_it32(i,t) = y_g3(i,t) - (c2) - z3(t,:,i)*pi_g3;
    end
end
mu_star_i32 = -s_u2^2/s_i32^2*sum(epsilon_it32,2); % n X 1 summation over t
s_star_i32 = sqrt(s_u2^2*s_v3^2/s_i32^2);

% part 1 of the loglikelihood
%P1 = -nn*(T-1)/2*log(s_v^2) ;

P1_g11 = -(T1-1)/2*log(s_v1^2);
P1_g21 = -(T2-1)/2*log(s_v2^2);
P1_g31 = -(T3-1)/2*log(s_v3^2);

P1_g12 = -(T1-1)/2*log(s_v1^2);
P1_g22 = -(T2-1)/2*log(s_v2^2);
P1_g32 = -(T3-1)/2*log(s_v3^2);

% part 2 
%P2 = -(nn/2)*log(s_i^2);

P2_g11 = -(1/2)*log(s_i11^2);
P2_g21 = -(1/2)*log(s_i21^2);
P2_g31 = -(1/2)*log(s_i31^2);

P2_g12 = -(1/2)*log(s_i12^2);
P2_g22 = -(1/2)*log(s_i22^2);
P2_g32 = -(1/2)*log(s_i32^2);

% part 3 
%P3 = sum(log(1-normcdf(-mu_star_i/s_star_i)));

P3_g11 = (log(1-normcdf(-mu_star_i11/s_star_i11)));
P3_g21 = (log(1-normcdf(-mu_star_i21/s_star_i21)));
P3_g31 = (log(1-normcdf(-mu_star_i31/s_star_i31)));

P3_g12 = (log(1-normcdf(-mu_star_i12/s_star_i12)));
P3_g22 = (log(1-normcdf(-mu_star_i22/s_star_i22)));
P3_g32 = (log(1-normcdf(-mu_star_i32/s_star_i32)));


% part 4 
%P4 = (1/2)*sum((mu_star_i/s_star_i).^2);

P4_g11 = (1/2)*(mu_star_i11/s_star_i11).^2;
P4_g21 = (1/2)*(mu_star_i21/s_star_i21).^2;
P4_g31 = (1/2)*(mu_star_i31/s_star_i31).^2;

P4_g12 = (1/2)*(mu_star_i12/s_star_i12).^2;
P4_g22 = (1/2)*(mu_star_i22/s_star_i22).^2;
P4_g32 = (1/2)*(mu_star_i32/s_star_i32).^2;

% part 5
%P5 = - 1/(2*s_v^2)*sum(sum(epsilon_it.^2));

P5_g11 = -1/(2*s_v1^2)*sum(epsilon_it11.^2,2);
P5_g21 = -1/(2*s_v2^2)*sum(epsilon_it21.^2,2);
P5_g31 = -1/(2*s_v3^2)*sum(epsilon_it31.^2,2);

P5_g12 = -1/(2*s_v1^2)*sum(epsilon_it12.^2,2);
P5_g22 = -1/(2*s_v2^2)*sum(epsilon_it22.^2,2);
P5_g32 = -1/(2*s_v3^2)*sum(epsilon_it32.^2,2);

% the minus log likelihood function
% L = -(P1 + P2 + P3 + P4 + P5);

L11 = -(P1_g11 + P2_g11 + P3_g11 + P4_g11 + P5_g11);
L12 = -(P1_g12 + P2_g12 + P3_g12 + P4_g12 + P5_g12);

L21 = -(P1_g21 + P2_g21 + P3_g21 + P4_g21 + P5_g21);
L22 = -(P1_g22 + P2_g22 + P3_g22 + P4_g22 + P5_g22);

L31 = -(P1_g31 + P2_g31 + P3_g31 + P4_g31 + P5_g31);
L32 = -(P1_g32 + P2_g32 + P3_g32 + P4_g32 + P5_g32);

L1 = sum(log(tau*exp(-L11) + (1-tau)*exp(-L12)));
L2 = sum(log(tau*exp(-L21) + (1-tau)*exp(-L22)));
L3 = sum(log(tau*exp(-L31) + (1-tau)*exp(-L32)));


Lprime = -(L1 + L2 + L3);
