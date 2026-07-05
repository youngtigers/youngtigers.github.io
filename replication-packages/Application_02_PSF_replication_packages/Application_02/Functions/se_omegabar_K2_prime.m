%%
%%  se_omegabar_prime =se_omegabar_prime_K2_prime(omegabar, y, z, sigvhat, pi) 
%%
%%  "loglik_K2_prime" computes negative of log-likelihood function
%%  used in step 4' of the estimation alogorithm for the case where
%%  c - u_i has a group structure, assuming there's 2 groups for
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
%%      se_omegabar_prime:  vector of se of omegabar prime evaluated at specfied values

function se_omegabar_prime = se_omegabar_K2_prime(omegabar, y, z1, z2, sigvhat_g1, sigvhat_g2, pi_g1, pi_g2, K_sieves ) 

% set numerical approx constant
h = 10^-5;

% Parameters: same ordering as loglik_K2_prime
c1        = omegabar(1);
sqrt_s_u1 = omegabar(2); % sqrt_s_u is the square root of sigma_u1; this ensures sigma_u1 is always positive
c2        = omegabar(3);
sqrt_s_u2 = omegabar(4); % sqrt_s_u is the square root of sigma_u2; this ensures sigma_u2 is always positive
tau       = omegabar(5);

% Parameters: Note the change of variable
sqrt_s_v1 = sqrt(sigvhat_g1); % sqrt_s_v1 is the square root of sigma_v (in std); note sigvhat = sqrt(sigv2hat_g1)
sqrt_s_v2 = sqrt(sigvhat_g2); % sqrt_s_v2 is the square root of sigma_v (in std); note sigvhat = sqrt(sigv2hat_g2)



% sort response and covariates into two groups based on HAC
y_g1  =  y(K_sieves==1,:);
[nn1,T1] = size(y_g1); % get the sample size
        
y_g2  =  y(K_sieves==2,:);
[nn2,T2] = size(y_g2); % get the sample size

N = nn1 + nn2;

% conversion
s_v1 = sqrt_s_v1^2;   % square it to make it sigma_v1
s_v2 = sqrt_s_v2^2;   % square it to make it sigma_v2
s_u1 = sqrt_s_u1^2;   % square it to make it sigma_u1
s_u2 = sqrt_s_u2^2;   % square it to make it sigma_u2

% set storage
epsilon_it11 = NaN(nn1,T1);
epsilon_it12 = NaN(nn1,T1);
epsilon_it21 = NaN(nn2,T2);
epsilon_it22 = NaN(nn2,T2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimate se(c1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for el = [1 0 -1]

% Group 1, a1, u1
s_i11 = sqrt(s_v1^2 + T1*s_u1^2);      
for i = 1:nn1
    for t = 1:T1
        epsilon_it11(i,t) = y_g1(i,t) - (c1 + h*el) - z1(t,:,i)*pi_g1;   % n X T the residuls
    end
end
mu_star_i11 = -s_u1^2/s_i11^2*sum(epsilon_it11,2); % n X 1 summation over t
s_star_i11 = sqrt(s_u1^2*s_v1^2/s_i11^2);

% Group 1, a2, u2
s_i12 = sqrt(s_v1^2 + T1*s_u2^2);
%epsilon_it12 = y_g1 - c2 - x1_g1*b1_g1 - x2_g1*b2_g1;         % n X T the residuls
for i = 1:nn1
    for t = 1:T1
        epsilon_it12(i,t) = y_g1(i,t) - c2 - z1(t,:,i)*pi_g1;
    end
end
mu_star_i12 = -s_u2^2/s_i12^2*sum(epsilon_it12,2); % n X 1 summation over t
s_star_i12 = sqrt(s_u2^2*s_v1^2/s_i12^2);

% Group 2,  a1, u1
s_i21 = sqrt(s_v2^2 + T2*s_u1^2);
%epsilon_it21 = y_g2 - a1 - x1_g2*b1_g2 - x2_g2*b2_g2;         % n X T the residuls
for i = 1:nn2
    for t = 1:T2
        epsilon_it21(i,t) = y_g2(i,t) - (c1 + h*el) - z2(t,:,i)*pi_g2;
    end
end
mu_star_i21 = -s_u1^2/s_i21^2*sum(epsilon_it21,2); % n X 1 summation over t
s_star_i21 = sqrt(s_u1^2*s_v2^2/s_i21^2);

% Group 2,  a2, u2
s_i22 = sqrt(s_v2^2 + T2*s_u2^2);
%epsilon_it22 = y_g2 - a2 - x1_g2*b1_g2 - x2_g2*b2_g2;         % n X T the residuls
for i = 1:nn2
    for t = 1:T2
        epsilon_it22(i,t) = y_g2(i,t) - c2 - z2(t,:,i)*pi_g2;
    end
end
mu_star_i22 = -s_u2^2/s_i22^2*sum(epsilon_it22,2); % n X 1 summation over t
s_star_i22 = sqrt(s_u2^2*s_v2^2/s_i22^2);

% part 1 of the loglikelihood
P1_g11 = -(T1-1)/2*log(s_v1^2);
P1_g12 = -(T1-1)/2*log(s_v1^2);

P1_g21 = -(T2-1)/2*log(s_v2^2);
P1_g22 = -(T2-1)/2*log(s_v2^2);

% part 2 
P2_g11 = -(1/2)*log(s_i11^2);
P2_g12 = -(1/2)*log(s_i12^2);

P2_g21 = -(1/2)*log(s_i21^2);
P2_g22 = -(1/2)*log(s_i22^2);

% part 3 
P3_g11 = log(1-normcdf(-mu_star_i11/s_star_i11));
P3_g12 = log(1-normcdf(-mu_star_i12/s_star_i12));

P3_g21 = log(1-normcdf(-mu_star_i21/s_star_i21));
P3_g22 = log(1-normcdf(-mu_star_i22/s_star_i22));

% part 4 
P4_g11 = (1/2)*(mu_star_i11/s_star_i11).^2;
P4_g12 = (1/2)*(mu_star_i12/s_star_i12).^2;

P4_g21 = (1/2)*(mu_star_i21/s_star_i21).^2;
P4_g22 = (1/2)*(mu_star_i22/s_star_i22).^2;

% part 5
%P5 = - 1/(2*s_v^2)*sum(sum(epsilon_it.^2));

P5_g11 = -1/(2*s_v1^2)*sum(epsilon_it11.^2,2);
P5_g12 = -1/(2*s_v1^2)*sum(epsilon_it12.^2,2);

P5_g21 = -1/(2*s_v2^2)*sum(epsilon_it21.^2,2);
P5_g22 = -1/(2*s_v2^2)*sum(epsilon_it22.^2,2);

% the minus log likelihood function
% L = -(P1 + P2 + P3 + P4 + P5);

L11 = (P1_g11 + P2_g11 + P3_g11 + P4_g11 + P5_g11); % group 1, mixture 1
L12 = (P1_g12 + P2_g12 + P3_g12 + P4_g12 + P5_g12); % group 1, mixture 2

L21 = (P1_g21 + P2_g21 + P3_g21 + P4_g21 + P5_g21); % group 2, mixture 1
L22 = (P1_g22 + P2_g22 + P3_g22 + P4_g22 + P5_g22); % group 2, mixture 2

L1 = sum(log(tau*exp(L11) + (1-tau)*exp(L12)));    % log-likelihood of group 1, summed over group 1 observations
L2 = sum(log(tau*exp(L21) + (1-tau)*exp(L22)));    % log-likelihood of group 2


if el == 1
    Lprime_c1(1) = (L1 + L2)/N;
elseif el == 0
    Lprime_c1(2) = (L1 + L2)/N;
elseif el == -1
    Lprime_c1(3) = (L1 + L2)/N;
end

end


II_c1 = -( (Lprime_c1(1)) - 2*(Lprime_c1(2)) + (Lprime_c1(3))) / (h^2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimate se(sigu1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for el = [1 0 -1]

% Group 1, a1, u1
s_i11 = sqrt(s_v1^2 + T1*(s_u1+ h*el)^2);      
for i = 1:nn1
    for t = 1:T1
        epsilon_it11(i,t) = y_g1(i,t) - c1 - z1(t,:,i)*pi_g1;   % n X T the residuls
    end
end
mu_star_i11 = -(s_u1+ h*el)^2/s_i11^2*sum(epsilon_it11,2); % n X 1 summation over t
s_star_i11 = sqrt((s_u1+ h*el)^2*s_v1^2/s_i11^2);

% Group 1, a2, u2
s_i12 = sqrt(s_v1^2 + T1*s_u2^2);
%epsilon_it12 = y_g1 - c2 - x1_g1*b1_g1 - x2_g1*b2_g1;         % n X T the residuls
for i = 1:nn1
    for t = 1:T1
        epsilon_it12(i,t) = y_g1(i,t) - c2 - z1(t,:,i)*pi_g1;
    end
end
mu_star_i12 = -s_u2^2/s_i12^2*sum(epsilon_it12,2); % n X 1 summation over t
s_star_i12 = sqrt(s_u2^2*s_v1^2/s_i12^2);

% Group 2,  a1, u1
s_i21 = sqrt(s_v2^2 + T2*(s_u1+ h*el)^2);
%epsilon_it21 = y_g2 - a1 - x1_g2*b1_g2 - x2_g2*b2_g2;         % n X T the residuls
for i = 1:nn2
    for t = 1:T2
        epsilon_it21(i,t) = y_g2(i,t) - c1 - z2(t,:,i)*pi_g2;
    end
end
mu_star_i21 = -(s_u1+ h*el)^2/s_i21^2*sum(epsilon_it21,2); % n X 1 summation over t
s_star_i21 = sqrt((s_u1+ h*el)^2*s_v2^2/s_i21^2);

% Group 2,  a2, u2
s_i22 = sqrt(s_v2^2 + T2*s_u2^2);
%epsilon_it22 = y_g2 - a2 - x1_g2*b1_g2 - x2_g2*b2_g2;         % n X T the residuls
for i = 1:nn2
    for t = 1:T2
        epsilon_it22(i,t) = y_g2(i,t) - c2 - z2(t,:,i)*pi_g2;
    end
end
mu_star_i22 = -s_u2^2/s_i22^2*sum(epsilon_it22,2); % n X 1 summation over t
s_star_i22 = sqrt(s_u2^2*s_v2^2/s_i22^2);

% part 1 of the loglikelihood
P1_g11 = -(T1-1)/2*log(s_v1^2);
P1_g12 = -(T1-1)/2*log(s_v1^2);

P1_g21 = -(T2-1)/2*log(s_v2^2);
P1_g22 = -(T2-1)/2*log(s_v2^2);

% part 2 
P2_g11 = -(1/2)*log(s_i11^2);
P2_g12 = -(1/2)*log(s_i12^2);

P2_g21 = -(1/2)*log(s_i21^2);
P2_g22 = -(1/2)*log(s_i22^2);

% part 3 
P3_g11 = log(1-normcdf(-mu_star_i11/s_star_i11));
P3_g12 = log(1-normcdf(-mu_star_i12/s_star_i12));

P3_g21 = log(1-normcdf(-mu_star_i21/s_star_i21));
P3_g22 = log(1-normcdf(-mu_star_i22/s_star_i22));

% part 4 
P4_g11 = (1/2)*(mu_star_i11/s_star_i11).^2;
P4_g12 = (1/2)*(mu_star_i12/s_star_i12).^2;

P4_g21 = (1/2)*(mu_star_i21/s_star_i21).^2;
P4_g22 = (1/2)*(mu_star_i22/s_star_i22).^2;

% part 5
%P5 = - 1/(2*s_v^2)*sum(sum(epsilon_it.^2));

P5_g11 = -1/(2*s_v1^2)*sum(epsilon_it11.^2,2);
P5_g12 = -1/(2*s_v1^2)*sum(epsilon_it12.^2,2);

P5_g21 = -1/(2*s_v2^2)*sum(epsilon_it21.^2,2);
P5_g22 = -1/(2*s_v2^2)*sum(epsilon_it22.^2,2);

% the minus log likelihood function
% L = -(P1 + P2 + P3 + P4 + P5);

L11 = (P1_g11 + P2_g11 + P3_g11 + P4_g11 + P5_g11); % group 1, mixture 1
L12 = (P1_g12 + P2_g12 + P3_g12 + P4_g12 + P5_g12); % group 1, mixture 2

L21 = (P1_g21 + P2_g21 + P3_g21 + P4_g21 + P5_g21); % group 2, mixture 1
L22 = (P1_g22 + P2_g22 + P3_g22 + P4_g22 + P5_g22); % group 2, mixture 2

L1 = sum(log(tau*exp(L11) + (1-tau)*exp(L12)));    % log-likelihood of group 1, summed over group 1 observations
L2 = sum(log(tau*exp(L21) + (1-tau)*exp(L22)));    % log-likelihood of group 2


if el == 1
    Lprime_sigu1(1) = (L1 + L2)/N;
elseif el == 0
    Lprime_sigu1(2) = (L1 + L2)/N;
elseif el == -1
    Lprime_sigu1(3) = (L1 + L2)/N;
end

end


II_sigu1 = -( (Lprime_sigu1(1)) - 2*(Lprime_sigu1(2)) + (Lprime_sigu1(3)) ) / (h^2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimate se(c2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for el = [1 0 -1]

% Group 1, a1, u1
s_i11 = sqrt(s_v1^2 + T1*s_u1^2);      
for i = 1:nn1
    for t = 1:T1
        epsilon_it11(i,t) = y_g1(i,t) - c1 - z1(t,:,i)*pi_g1;   % n X T the residuls
    end
end
mu_star_i11 = -s_u1^2/s_i11^2*sum(epsilon_it11,2); % n X 1 summation over t
s_star_i11 = sqrt(s_u1^2*s_v1^2/s_i11^2);

% Group 1, a2, u2
s_i12 = sqrt(s_v1^2 + T1*s_u2^2);
%epsilon_it12 = y_g1 - c2 - x1_g1*b1_g1 - x2_g1*b2_g1;         % n X T the residuls
for i = 1:nn1
    for t = 1:T1
        epsilon_it12(i,t) = y_g1(i,t) - (c2+ h*el) - z1(t,:,i)*pi_g1;
    end
end
mu_star_i12 = -s_u2^2/s_i12^2*sum(epsilon_it12,2); % n X 1 summation over t
s_star_i12 = sqrt(s_u2^2*s_v1^2/s_i12^2);

% Group 2,  a1, u1
s_i21 = sqrt(s_v2^2 + T2*s_u1^2);
%epsilon_it21 = y_g2 - a1 - x1_g2*b1_g2 - x2_g2*b2_g2;         % n X T the residuls
for i = 1:nn2
    for t = 1:T2
        epsilon_it21(i,t) = y_g2(i,t) - c1 - z2(t,:,i)*pi_g2;
    end
end
mu_star_i21 = -s_u1^2/s_i21^2*sum(epsilon_it21,2); % n X 1 summation over t
s_star_i21 = sqrt(s_u1^2*s_v2^2/s_i21^2);

% Group 2,  a2, u2
s_i22 = sqrt(s_v2^2 + T2*s_u2^2);
%epsilon_it22 = y_g2 - a2 - x1_g2*b1_g2 - x2_g2*b2_g2;         % n X T the residuls
for i = 1:nn2
    for t = 1:T2
        epsilon_it22(i,t) = y_g2(i,t) - (c2+ h*el) - z2(t,:,i)*pi_g2;
    end
end
mu_star_i22 = -s_u2^2/s_i22^2*sum(epsilon_it22,2); % n X 1 summation over t
s_star_i22 = sqrt(s_u2^2*s_v2^2/s_i22^2);

% part 1 of the loglikelihood
P1_g11 = -(T1-1)/2*log(s_v1^2);
P1_g12 = -(T1-1)/2*log(s_v1^2);

P1_g21 = -(T2-1)/2*log(s_v2^2);
P1_g22 = -(T2-1)/2*log(s_v2^2);

% part 2 
P2_g11 = -(1/2)*log(s_i11^2);
P2_g12 = -(1/2)*log(s_i12^2);

P2_g21 = -(1/2)*log(s_i21^2);
P2_g22 = -(1/2)*log(s_i22^2);

% part 3 
P3_g11 = log(1-normcdf(-mu_star_i11/s_star_i11));
P3_g12 = log(1-normcdf(-mu_star_i12/s_star_i12));

P3_g21 = log(1-normcdf(-mu_star_i21/s_star_i21));
P3_g22 = log(1-normcdf(-mu_star_i22/s_star_i22));

% part 4 
P4_g11 = (1/2)*(mu_star_i11/s_star_i11).^2;
P4_g12 = (1/2)*(mu_star_i12/s_star_i12).^2;

P4_g21 = (1/2)*(mu_star_i21/s_star_i21).^2;
P4_g22 = (1/2)*(mu_star_i22/s_star_i22).^2;

% part 5
%P5 = - 1/(2*s_v^2)*sum(sum(epsilon_it.^2));

P5_g11 = -1/(2*s_v1^2)*sum(epsilon_it11.^2,2);
P5_g12 = -1/(2*s_v1^2)*sum(epsilon_it12.^2,2);

P5_g21 = -1/(2*s_v2^2)*sum(epsilon_it21.^2,2);
P5_g22 = -1/(2*s_v2^2)*sum(epsilon_it22.^2,2);

% the minus log likelihood function
% L = -(P1 + P2 + P3 + P4 + P5);

L11 = (P1_g11 + P2_g11 + P3_g11 + P4_g11 + P5_g11); % group 1, mixture 1
L12 = (P1_g12 + P2_g12 + P3_g12 + P4_g12 + P5_g12); % group 1, mixture 2

L21 = (P1_g21 + P2_g21 + P3_g21 + P4_g21 + P5_g21); % group 2, mixture 1
L22 = (P1_g22 + P2_g22 + P3_g22 + P4_g22 + P5_g22); % group 2, mixture 2

L1 = sum(log(tau*exp(L11) + (1-tau)*exp(L12)));    % log-likelihood of group 1, summed over group 1 observations
L2 = sum(log(tau*exp(L21) + (1-tau)*exp(L22)));    % log-likelihood of group 2


if el == 1
    Lprime_c2(1) = (L1 + L2)/N;
elseif el == 0
    Lprime_c2(2) = (L1 + L2)/N;
elseif el == -1
    Lprime_c2(3) = (L1 + L2)/N;
end

end


II_c2 = -( (Lprime_c2(1)) - 2*(Lprime_c2(2)) + (Lprime_c2(3))) / (h^2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimate se(sigu2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for el = [1 0 -1]

% Group 1, a1, u1
s_i11 = sqrt(s_v1^2 + T1*s_u1^2);      
for i = 1:nn1
    for t = 1:T1
        epsilon_it11(i,t) = y_g1(i,t) - c1 - z1(t,:,i)*pi_g1;   % n X T the residuls
    end
end
mu_star_i11 = -s_u1^2/s_i11^2*sum(epsilon_it11,2); % n X 1 summation over t
s_star_i11 = sqrt(s_u1^2*s_v1^2/s_i11^2);

% Group 1, a2, u2
s_i12 = sqrt(s_v1^2 + T1*(s_u2+ h*el)^2);
%epsilon_it12 = y_g1 - c2 - x1_g1*b1_g1 - x2_g1*b2_g1;         % n X T the residuls
for i = 1:nn1
    for t = 1:T1
        epsilon_it12(i,t) = y_g1(i,t) - c2 - z1(t,:,i)*pi_g1;
    end
end
mu_star_i12 = -(s_u2+ h*el)^2/s_i12^2*sum(epsilon_it12,2); % n X 1 summation over t
s_star_i12 = sqrt((s_u2+ h*el)^2*s_v1^2/s_i12^2);

% Group 2,  a1, u1
s_i21 = sqrt(s_v2^2 + T2*s_u1^2);
%epsilon_it21 = y_g2 - a1 - x1_g2*b1_g2 - x2_g2*b2_g2;         % n X T the residuls
for i = 1:nn2
    for t = 1:T2
        epsilon_it21(i,t) = y_g2(i,t) - c1 - z2(t,:,i)*pi_g2;
    end
end
mu_star_i21 = -s_u1^2/s_i21^2*sum(epsilon_it21,2); % n X 1 summation over t
s_star_i21 = sqrt(s_u1^2*s_v2^2/s_i21^2);

% Group 2,  a2, u2
s_i22 = sqrt(s_v2^2 + T2*(s_u2+ h*el)^2);
%epsilon_it22 = y_g2 - a2 - x1_g2*b1_g2 - x2_g2*b2_g2;         % n X T the residuls
for i = 1:nn2
    for t = 1:T2
        epsilon_it22(i,t) = y_g2(i,t) - c2 - z2(t,:,i)*pi_g2;
    end
end
mu_star_i22 = -(s_u2+ h*el)^2/s_i22^2*sum(epsilon_it22,2); % n X 1 summation over t
s_star_i22 = sqrt((s_u2+ h*el)^2*s_v2^2/s_i22^2);

% part 1 of the loglikelihood
P1_g11 = -(T1-1)/2*log(s_v1^2);
P1_g12 = -(T1-1)/2*log(s_v1^2);

P1_g21 = -(T2-1)/2*log(s_v2^2);
P1_g22 = -(T2-1)/2*log(s_v2^2);

% part 2 
P2_g11 = -(1/2)*log(s_i11^2);
P2_g12 = -(1/2)*log(s_i12^2);

P2_g21 = -(1/2)*log(s_i21^2);
P2_g22 = -(1/2)*log(s_i22^2);

% part 3 
P3_g11 = log(1-normcdf(-mu_star_i11/s_star_i11));
P3_g12 = log(1-normcdf(-mu_star_i12/s_star_i12));

P3_g21 = log(1-normcdf(-mu_star_i21/s_star_i21));
P3_g22 = log(1-normcdf(-mu_star_i22/s_star_i22));

% part 4 
P4_g11 = (1/2)*(mu_star_i11/s_star_i11).^2;
P4_g12 = (1/2)*(mu_star_i12/s_star_i12).^2;

P4_g21 = (1/2)*(mu_star_i21/s_star_i21).^2;
P4_g22 = (1/2)*(mu_star_i22/s_star_i22).^2;

% part 5
%P5 = - 1/(2*s_v^2)*sum(sum(epsilon_it.^2));

P5_g11 = -1/(2*s_v1^2)*sum(epsilon_it11.^2,2);
P5_g12 = -1/(2*s_v1^2)*sum(epsilon_it12.^2,2);

P5_g21 = -1/(2*s_v2^2)*sum(epsilon_it21.^2,2);
P5_g22 = -1/(2*s_v2^2)*sum(epsilon_it22.^2,2);

% the minus log likelihood function
% L = -(P1 + P2 + P3 + P4 + P5);

L11 = (P1_g11 + P2_g11 + P3_g11 + P4_g11 + P5_g11); % group 1, mixture 1
L12 = (P1_g12 + P2_g12 + P3_g12 + P4_g12 + P5_g12); % group 1, mixture 2

L21 = (P1_g21 + P2_g21 + P3_g21 + P4_g21 + P5_g21); % group 2, mixture 1
L22 = (P1_g22 + P2_g22 + P3_g22 + P4_g22 + P5_g22); % group 2, mixture 2

L1 = sum(log(tau*exp(L11) + (1-tau)*exp(L12)));    % log-likelihood of group 1, summed over group 1 observations
L2 = sum(log(tau*exp(L21) + (1-tau)*exp(L22)));    % log-likelihood of group 2


if el == 1
    Lprime_sigu2(1) = (L1 + L2)/N;
elseif el == 0
    Lprime_sigu2(2) = (L1 + L2)/N;
elseif el == -1
    Lprime_sigu2(3) = (L1 + L2)/N;
end

end


II_sigu2 = -( (Lprime_sigu2(1)) - 2*(Lprime_sigu2(2)) + (Lprime_sigu2(3)) ) / (h^2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimate se(tau)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for el = [1 0 -1]

% Group 1, a1, u1
s_i11 = sqrt(s_v1^2 + T1*s_u1^2);      
for i = 1:nn1
    for t = 1:T1
        epsilon_it11(i,t) = y_g1(i,t) - c1 - z1(t,:,i)*pi_g1;   % n X T the residuls
    end
end
mu_star_i11 = -s_u1^2/s_i11^2*sum(epsilon_it11,2); % n X 1 summation over t
s_star_i11 = sqrt(s_u1^2*s_v1^2/s_i11^2);

% Group 1, a2, u2
s_i12 = sqrt(s_v1^2 + T1*s_u2^2);
%epsilon_it12 = y_g1 - c2 - x1_g1*b1_g1 - x2_g1*b2_g1;         % n X T the residuls
for i = 1:nn1
    for t = 1:T1
        epsilon_it12(i,t) = y_g1(i,t) - c2 - z1(t,:,i)*pi_g1;
    end
end
mu_star_i12 = -s_u2^2/s_i12^2*sum(epsilon_it12,2); % n X 1 summation over t
s_star_i12 = sqrt(s_u2^2*s_v1^2/s_i12^2);

% Group 2,  a1, u1
s_i21 = sqrt(s_v2^2 + T2*s_u1^2);
%epsilon_it21 = y_g2 - a1 - x1_g2*b1_g2 - x2_g2*b2_g2;         % n X T the residuls
for i = 1:nn2
    for t = 1:T2
        epsilon_it21(i,t) = y_g2(i,t) - c1 - z2(t,:,i)*pi_g2;
    end
end
mu_star_i21 = -s_u1^2/s_i21^2*sum(epsilon_it21,2); % n X 1 summation over t
s_star_i21 = sqrt(s_u1^2*s_v2^2/s_i21^2);

% Group 2,  a2, u2
s_i22 = sqrt(s_v2^2 + T2*s_u2^2);
%epsilon_it22 = y_g2 - a2 - x1_g2*b1_g2 - x2_g2*b2_g2;         % n X T the residuls
for i = 1:nn2
    for t = 1:T2
        epsilon_it22(i,t) = y_g2(i,t) - c2 - z2(t,:,i)*pi_g2;
    end
end
mu_star_i22 = -s_u2^2/s_i22^2*sum(epsilon_it22,2); % n X 1 summation over t
s_star_i22 = sqrt(s_u2^2*s_v2^2/s_i22^2);

% part 1 of the loglikelihood
P1_g11 = -(T1-1)/2*log(s_v1^2);
P1_g12 = -(T1-1)/2*log(s_v1^2);

P1_g21 = -(T2-1)/2*log(s_v2^2);
P1_g22 = -(T2-1)/2*log(s_v2^2);

% part 2 
P2_g11 = -(1/2)*log(s_i11^2);
P2_g12 = -(1/2)*log(s_i12^2);

P2_g21 = -(1/2)*log(s_i21^2);
P2_g22 = -(1/2)*log(s_i22^2);

% part 3 
P3_g11 = log(1-normcdf(-mu_star_i11/s_star_i11));
P3_g12 = log(1-normcdf(-mu_star_i12/s_star_i12));

P3_g21 = log(1-normcdf(-mu_star_i21/s_star_i21));
P3_g22 = log(1-normcdf(-mu_star_i22/s_star_i22));

% part 4 
P4_g11 = (1/2)*(mu_star_i11/s_star_i11).^2;
P4_g12 = (1/2)*(mu_star_i12/s_star_i12).^2;

P4_g21 = (1/2)*(mu_star_i21/s_star_i21).^2;
P4_g22 = (1/2)*(mu_star_i22/s_star_i22).^2;

% part 5
%P5 = - 1/(2*s_v^2)*sum(sum(epsilon_it.^2));

P5_g11 = -1/(2*s_v1^2)*sum(epsilon_it11.^2,2);
P5_g12 = -1/(2*s_v1^2)*sum(epsilon_it12.^2,2);

P5_g21 = -1/(2*s_v2^2)*sum(epsilon_it21.^2,2);
P5_g22 = -1/(2*s_v2^2)*sum(epsilon_it22.^2,2);

% the minus log likelihood function
% L = -(P1 + P2 + P3 + P4 + P5);

L11 = (P1_g11 + P2_g11 + P3_g11 + P4_g11 + P5_g11); % group 1, mixture 1
L12 = (P1_g12 + P2_g12 + P3_g12 + P4_g12 + P5_g12); % group 1, mixture 2

L21 = (P1_g21 + P2_g21 + P3_g21 + P4_g21 + P5_g21); % group 2, mixture 1
L22 = (P1_g22 + P2_g22 + P3_g22 + P4_g22 + P5_g22); % group 2, mixture 2

L1 = sum(log( (tau + h*el)*exp(L11) + (1-(tau+ h*el))*exp(L12)));    % log-likelihood of group 1, summed over group 1 observations
L2 = sum(log( (tau + h*el)*exp(L21) + (1-(tau+ h*el))*exp(L22)));    % log-likelihood of group 2

%L1 = sum(log( (tau + h*el)*exp(L12) + (1-(tau+ h*el))*exp(L11)));    % log-likelihood of group 1, summed over group 1 observations
%L2 = sum(log( (tau + h*el)*exp(L22) + (1-(tau+ h*el))*exp(L21)));    % log-likelihood of group 2


if el == 1
    Lprime_tau(1) = ( L1 + L2)/N;
elseif el == 0
    Lprime_tau(2) = ( L1 + L2)/N;
elseif el == -1
    Lprime_tau(3) = ( L1 + L2)/N;
end

end


II_tau = -( (Lprime_tau(1)) - 2*(Lprime_tau(2)) + (Lprime_tau(3))) / (h^2);


% Output order [tau, c1, sigu1, c2, sigu2] matching the display remap:
% II_c2 = curvature w.r.t. pos3 = paper's c1
% II_sigu2 = curvature w.r.t. pos4 = paper's sigu1
% II_c1 = curvature w.r.t. pos1 = paper's c2
% II_sigu1 = curvature w.r.t. pos2 = paper's sigu2
se_omegabar_prime = 1 ./ sqrt( N * [II_tau, II_c2, II_sigu2, II_c1, II_sigu1] );
