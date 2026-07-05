%%
%%  L = loglik_K1(omegabar, y, z, sigvhat, pi) 
%%
%%  "loglik_K1" computes negative of log-likelihood function
%%  used in step 4 of the estimation alogorithm for the case where
%%  c - u_i has NO group structure, assuming there's 1 group for
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


function L = loglik_K1(omegabar, y, z, sigvhat, pi) 

% get the sample size
[nn,T] = size(y); 

% initialize
epsilon_it = NaN(nn,T);

% Parameters
c = omegabar(1);
sqrt_s_u = omegabar(2); % sqrt_s_u is the square root of sigma_u; this ensures sigma_u is always positive
s_u = sqrt_s_u^2;   % square it to make it sigma_u

% estimated parameters
sqrt_s_v = sqrt(sigvhat);
s_v = sqrt_s_v^2;

% No grouping
s_i = sqrt(s_v^2 + T*s_u^2);
for i = 1:nn
    for t = 1:T
        epsilon_it(i,t) = y(i,t) - c - z(t,:,i)*pi;
    end
end
mu_star_i = -s_u^2/s_i^2*sum(epsilon_it,2); % n X 1 summation over t
s_star_i = sqrt(s_u^2*s_v^2/s_i^2);

% part 1 of the loglikelihood
P1 = -nn*(T-1)/2*log(s_v^2) ;

% part 2 
P2 = -(nn/2)*log(s_i^2);

% part 3 
P3 = sum(log(1-normcdf(-mu_star_i/s_star_i)));

% part 4 
P4 = (1/2)*sum((mu_star_i/s_star_i).^2);

% part 5
P5 = - 1/(2*s_v^2)*sum(sum(epsilon_it.^2));

% the minus log likelihood function
L = -(P1 + P2 + P3 + P4 + P5);


end