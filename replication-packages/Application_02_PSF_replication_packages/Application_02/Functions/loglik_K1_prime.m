%%
%%  Lprime = loglik_K1_prime(omegabar, y, z, sigvhat, pi) 
%%
%%  "loglik_K1_prime" computes negative of log-likelihood function
%%  used in step 4' of the estimation alogorithm for the case where
%%  c - u_i has a group structure,assuming there's 1 group for
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

function Lprime = loglik_K1_prime(omegabar, y, z, sigvhat, pi) 

% get the sample size
[nn,T] = size(y); 

% initialize
epsilon_it1 = NaN(nn,T);
epsilon_it2 = NaN(nn,T);

% Parameters
c1        = omegabar(1);
sqrt_s_u1 = omegabar(2); % sqrt_s_u is the square root of sigma_u1; this ensures sigma_u1 is always positive
c2        = omegabar(3);
sqrt_s_u2 = omegabar(4); % sqrt_s_u is the square root of sigma_u2; this ensures sigma_u2 is always positive
tau       = omegabar(5);

% estimated parameters
sqrt_s_v = sqrt(sigvhat);
s_v = sqrt_s_v^2;

s_u1 = sqrt_s_u1^2;      % square it to make it sigma_u1 (std)
s_u2 = sqrt_s_u2^2;      % square it to make it sigma_u2 (std)

% No grouping, a1, u1
s_i1 = sqrt(s_v^2 + T*s_u1^2);
for i = 1:nn
    for t = 1:T
        epsilon_it1(i,t) = y(i,t) - c1 - z(t,:,i)*pi;
    end
end
mu_star_i1 = -s_u1^2/s_i1^2*sum(epsilon_it1,2); % n X 1 summation over t
s_star_i1 = sqrt(s_u1^2*s_v^2/s_i1^2);

% No grouping, a2, u2
s_i2 = sqrt(s_v^2 + T*s_u2^2);
for i = 1:nn
    for t = 1:T
        epsilon_it2(i,t) = y(i,t) - c2 - z(t,:,i)*pi;
    end
end
mu_star_i2 = -s_u2^2/s_i2^2*sum(epsilon_it2,2); % n X 1 summation over t
s_star_i2 = sqrt(s_u2^2*s_v^2/s_i2^2);

% part 1 of the loglikelihood
P1 = -(T-1)/2*log(s_v^2) ;

% part 2 
P2_1 = -(1/2)*log(s_i1^2);
P2_2 = -(1/2)*log(s_i2^2);

% part 3 
P3_1 = log(1-normcdf(-mu_star_i1./s_star_i1));  
P3_2 = log(1-normcdf(-mu_star_i2./s_star_i2)); % returns -Inf: normcdf = 1, (./.)= +large

% part 4 
P4_1 = (1/2)*(mu_star_i1/s_star_i1).^2;
P4_2 = (1/2)*(mu_star_i2/s_star_i2).^2;

% part 5
P5_1 = -1/(2*s_v^2)*sum(epsilon_it1.^2,2);
P5_2 = -1/(2*s_v^2)*sum(epsilon_it2.^2,2);

% the minus log likelihood function
L1 = -(P1 + P2_1 + P3_1 + P4_1 + P5_1);
L2 = -(P1 + P2_2 + P3_2 + P4_2 + P5_2);

% the minus log likelihood function
Lprime = -sum(log(tau*(exp(-L1)) + (1-tau)*(exp(-L2))));

