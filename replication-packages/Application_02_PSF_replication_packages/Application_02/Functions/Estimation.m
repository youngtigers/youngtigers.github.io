function out = Estimation(Params, y, x1, x2, x3, x4, x5, x1_cl, x2_cl, x3_cl, x4_cl, x5_cl)
%% Estimation  Core estimation engine for the PSF banking application.
%
%  out = Estimation(Params, y, x1..x5, x1_cl..x5_cl)
%
%  Runs the full estimation algorithm (Step 1 OLS sieve estimation, the
%  agglomerative clustering / information-criterion group selection, and the
%  Step 4/4' ML estimation of c, sigma_u and the mixture weight tau). This is
%  the code that previously lived inline in main_application.m; it is factored
%  out so that main_application.m is a thin driver.
%
%  Every variable the algorithm creates is returned, bundled, in the struct
%  OUT. The driver unpacks OUT back into the base workspace so the downstream
%  scripts (getFrontiers, Asymptotics, Frontiers_SE, Figures_Application,
%  sigv_omegabar_estimates) see exactly the variables they did before.
%
%  Inputs:
%    Params              settings/parameter struct (see main_application.m)
%    y                   N x T response (log cost)
%    x1..x5              N x T raw (log) regressors
%    x1_cl..x5_cl        N x T cleaned (demeaned/standardized) regressors
%
%  Output:
%    out                 struct holding every variable created below.

% Unpack the scalars the algorithm uses (same as main_application.m lines 60-66)
T         = Params.T;
N         = Params.N;
p         = Params.p;
m         = Params.m;
K_group   = Params.K_group;
lambda_s3 = Params.lambda_s3;
mLbar_S3  = Params.mLbar_S3;

%% Algorithm
fprintf('\n %%-- running algorithm --%% \n');
% Set storage
SetStorage

% ------------------------------- %
% Step 1: separate OLS estimation %
% ------------------------------- % 

% Cosine Basis functions {B_j}, column vectors
[B0, B1, B2, B3, B4, B5, B6, B7, B8, B9, B10, B11] = GenBasis(Params);
   
% construct Bm := \mathbb{B}^{m} 
if m >= 2
    Bm = GenBm(B1,B2,B3,B4,B5,B6,B7,B8,B9,B10,B11,Params); 
end

% construct xBm := x_it \kron \mathbb{B}^{m}
xBm = GenxBm(x1_cl,x2_cl,x3_cl,x4_cl,x5_cl,B1,B2,B3,B4,B5,B6,B7,B8,B9,B10,B11,Params);

% initialize
Zim = NaN(T,m*(p+1),N); 
pihat = NaN(m*(p+1),N);
yhat = NaN(N,T);
res = NaN(N,T);     % residual of the model in step 1


% construct Zim := Z_{im} and estimate pihat := \hat{\tilde{\pi}}_{i}
% (although pihat should really be pihat := \hat{\tilde{\pi}}_{i}\setminus{first element})
for i = 1:N
    if m == 1
        %Zim(:,:,i) = [ones(T,1) x1(i,:)', x2(i,:)' x3(i,:)' x4(i,:)'];
        Zim(:,:,i) = [ones(T,1) xBm(1:T,:,i)];
    else
        Zim(:,:,i) = [ones(T,1) Bm(1:T,:) xBm(1:T,:,i)];
    end
    pihat(:,i) = pinv(Zim(:,:,i)'*Zim(:,:,i))*Zim(:,:,i)'*y(i,:)';
end 

% calculate residual (better way)
for i = 1:N
    yhat(i,:) = Zim(:,:,i)*pihat(:,i);
    res(i,:) = ( y(i,:)' - Zim(:,:,i)*pihat(:,i) );
end 

% estimate sigvhat2 := \sigma_{vi}^{2} 
sigvhat2 = sum(res.^2,2) ./ (T-1);

% estimate R^2::
%Rsquare = 1 - sigvhat2'./var(y');
%fprintf('\n R^2 quantiles  (0.2, 0.5, 0.8): %.4f   %.4f   %.4f \n', ...
%    quantile(Rsquare, 0.2), quantile(Rsquare, 0.5), quantile(Rsquare, 0.8));

% define vthetahat := [pihat', \sigvhat2]' (exclude first element)
vthetahat = [pihat(2:end,:)', sigvhat2]';  % remove first element of pihat
    

% ------------------------------------------------- %
% Step 2: agglomerative hierarchical cluster tree Z %
% ------------------------------------------------- %

% Z encodes a tree containing hierarchical clusters of the rows of the input
Z_sieves = linkage(vthetahat(1:end,:)','ward','euclidean'); % group based on vthetahat    
   

for k = 1 : K_group

    % cluster parameters into max of k groups, each row of k is # group firm i belongs to
    K_sieves = cluster(Z_sieves,'maxclust',k);      % max number of cluster/group is ,K_group
    store_group_sieves(:,1,k) = K_sieves;
    
    % -------------------------------------- %
    % step 3: Post classification estimation %
    % -------------------------------------- %

    if k == 1

        % construct response and regressors 
        if isnan(mLbar_S3)
            mLbar_K1 = floor((N*T)^0.2); % set m lower bar
        else
            mLbar_K1 = mLbar_S3;         % set m lower bar
        end
        
        ZimLbar_K1 = inputs_g0(y,x1,x2,x3,x4,x5,B1, B2, B3, B4, B5, B6,B7,B8,B9,B10,B11, mLbar_K1, Params);
        zdot    = ZimLbar_K1 - mean(ZimLbar_K1,1);
        ydot    = y - mean(y,2);

        % stack response and regressors into matrix form
        % permute: convert dim(zdot) = T by mLbar * (p+1)-1 by N to T by N by mLbar * (p+1)-1
        % reshape: convert 3D object zdot into 2D object of size N*T by mLbar * (p+1)-1
        zdotK1 = reshape(permute(zdot,[1 3 2]), N*T, []);   
        % reshape: stack ydot matrix of size N by T to N*T by 1 vector
        ydotK1 = reshape(ydot', N*T, 1);

        % find pi via OLS
        pi_K1 = pinv(zdotK1'*zdotK1)*zdotK1'*ydotK1;

        % residual sum of squares
        RSS1 = (ydotK1 - zdotK1*pi_K1)'*(ydotK1 - zdotK1*pi_K1);
       
        % estimate variance of vit
        sigv2hat_K1 = (1/(N*(T-1)))*RSS1;
        
        % estimate standard error of sigv2hat
        vit = ydotK1 - zdotK1*pi_K1;
        vit2bar = sum(vit.^2)/(N*T);
        varhatvit2 = (1/(N*(T-1)))*sum( (vit.^2 - vit2bar).^2 );
        SE_sigvhat_K1 = ( sqrt(varhatvit2)/sqrt(N*T) )/(2*sqrt(sigv2hat_K1));

        % information criterion
        % No lambda
        IC_k1_L0 = N*T*log(sqrt(sigv2hat_K1)) + RSS1./sigv2hat_K1;

        IC_k1_L1 = IC_k1_L0 + lambda_s3*k;    % c_lambda=1   
        
    elseif k == 2
        
        %K_sieves = abs(K_sieves-3);
        
        % construct group 1 response and regressors 
        y_g1  =  y(K_sieves==1,:);
        [N_g1,T1] = size(y_g1); 
        if isnan(mLbar_S3)
            mLbar_K2g1 = floor((N_g1*T)^0.2); % set m lower bar
        else
            mLbar_K2g1 = mLbar_S3;
        end
        ZimLbar_K2_g1 = inputs_g(y,x1,x2,x3,x4,x5,B1, B2, B3, B4, B5,B6, B7, B8, B9, B10, B11,K_sieves,1,mLbar_K2g1,Params);       
        ydot_g1 = y_g1 - mean(y_g1,2);
        zdot_g1 = ZimLbar_K2_g1 - mean(ZimLbar_K2_g1,1);
        
        % construct group 2 response and regressors
        y_g2  =  y(K_sieves==2,:);
        [N_g2,T2] = size(y_g2);
        if isnan(mLbar_S3)
            mLbar_K2g2 = floor((N_g2*T)^0.2);   % set m lower bar
        else
            mLbar_K2g2 = mLbar_S3;
        end
        
        ZimLbar_K2_g2 = inputs_g(y,x1,x2,x3,x4,x5,B1, B2, B3, B4, B5, B6, B7, B8, B9, B10, B11, K_sieves,2,mLbar_K2g2,Params);        
        ydot_g2 = y_g2 - mean(y_g2,2);
        zdot_g2 = ZimLbar_K2_g2 - mean(ZimLbar_K2_g2,1);

        % stack response and regressors into matrix form
        % permute: convert dim(zdot) = T by mLbar * (p+1)-1 by N to T by N by mLbar * (p+1)-1
        % reshape: convert 3D object zdot into 2D object of size N*T by mLbar * (p+1)-1
        zdot_K2_g1 = reshape(permute(zdot_g1,[1 3 2]), N_g1*T, []);
        % reshape: stack ydot matrix of size N by T to N*T by 1 vector
        ydot_K2_g1 = reshape(ydot_g1',N_g1*T,1);
        zdot_K2_g2 = reshape(permute(zdot_g2,[1 3 2]), N_g2*T, []);
        ydot_K2_g2 = reshape(ydot_g2',N_g2*T,1);

        % find pi = argmin of RSS        
        pi_K2_g1 = pinv(zdot_K2_g1'*zdot_K2_g1)*zdot_K2_g1'*ydot_K2_g1;
        pi_K2_g2 = pinv(zdot_K2_g2'*zdot_K2_g2)*zdot_K2_g2'*ydot_K2_g2;

        % residual sum of squares
        RSS1 = (ydot_K2_g1 - zdot_K2_g1*pi_K2_g1)'*(ydot_K2_g1 - zdot_K2_g1*pi_K2_g1);
        RSS2 = (ydot_K2_g2 - zdot_K2_g2*pi_K2_g2)'*(ydot_K2_g2 - zdot_K2_g2*pi_K2_g2);

        % estimate variance of vit
        sigv2hat_K2_g1 = ((1/(N_g1*(T-1)))*RSS1);
        sigv2hat_K2_g2 = ((1/(N_g2*(T-1)))*RSS2);
        
        % estimate standard error of sigv2hat
        vit_g1 = ydot_K2_g1 - zdot_K2_g1*pi_K2_g1;
        vit_g2 = ydot_K2_g2 - zdot_K2_g2*pi_K2_g2;

        vit2bar_g1 = sum(vit_g1.^2)/(N_g1*T);
        vit2bar_g2 = sum(vit_g2.^2)/(N_g2*T);

        varhatvit2_g1 = (1/(N_g1*(T-1)))*sum( (vit_g1.^2 - vit2bar_g1).^2 );
        varhatvit2_g2 = (1/(N_g2*(T-1)))*sum( (vit_g2.^2 - vit2bar_g2).^2 );

        SE_sigvhat_K2_g1 = ( sqrt(varhatvit2_g1)/sqrt(N_g1*T) )/(2*sqrt(sigv2hat_K2_g1));
        SE_sigvhat_K2_g2 = ( sqrt(varhatvit2_g2)/sqrt(N_g2*T) )/(2*sqrt(sigv2hat_K2_g2));

        % information criterion
        % No lambda
        IC_k2_L0 = N_g1*T*log(sqrt(sigv2hat_K2_g1)) + RSS1./sigv2hat_K2_g1 ...
                 + N_g2*T*log(sqrt(sigv2hat_K2_g2)) + RSS2./sigv2hat_K2_g2;

        IC_k2_L1 = IC_k2_L0 + lambda_s3*k;    % c_lambda=1
        
        fprintf('size of clusters if K*=2: %s\n', ...
            mat2str([size(ydot_g1,1),size(ydot_g2,1)]));

    elseif k == 3 
        % construct group 1 response and regressors 
        y_g1  =  y(K_sieves==1,:);
        [N_g1,T1] = size(y_g1);
        if isnan(mLbar_S3)
            mLbar_K3g1 = floor((N_g1*T)^0.2);   % set m lower bar
        else
            mLbar_K3g1 = mLbar_S3;
        end
        ZimLbar_K3_g1 = inputs_g(y,x1,x2,x3,x4,x5,B1, B2, B3, B4, B5, B6, B7, B8, B9, B10, B11, K_sieves,1,mLbar_K3g1,Params);        
        ydot_g1 = y_g1 - mean(y_g1,2);
        zdot_g1 = ZimLbar_K3_g1 - mean(ZimLbar_K3_g1,1);
        
        % construct group 2 response and regressors
        y_g2  =  y(K_sieves==2,:);
        [N_g2,T2] = size(y_g2);
        if isnan(mLbar_S3)
            mLbar_K3g2 = floor((N_g2*T)^0.2);   % set m lower bar
        else
        mLbar_K3g2 = mLbar_S3;
        end
        ZimLbar_K3_g2 = inputs_g(y,x1,x2,x3,x4,x5,B1, B2, B3, B4, B5, B6, B7, B8, B9, B10, B11, K_sieves,2,mLbar_K3g2,Params);        
        ydot_g2 = y_g2 - mean(y_g2,2);
        zdot_g2 = ZimLbar_K3_g2 - mean(ZimLbar_K3_g2,1);
        
        % construct group 3 response and regressors
        y_g3  =  y(K_sieves==3,:);
        [N_g3,T3] = size(y_g3);
        if isnan(mLbar_S3)
            mLbar_K3g3 = floor((N_g3*T)^0.2);   % set m lower bar
        else
            mLbar_K3g3 = mLbar_S3;
        end
        ZimLbar_K3_g3 = inputs_g(y,x1,x2,x3,x4,x5,B1, B2, B3, B4, B5, B6, B7, B8, B9, B10, B11, K_sieves,3,mLbar_K3g3,Params);        
        ydot_g3 = y_g3 - mean(y_g3,2);
        zdot_g3 = ZimLbar_K3_g3 - mean(ZimLbar_K3_g3,1);
        
        % stack response and regressors into matrix form
        zdot_K3_g1 = reshape(permute(zdot_g1,[1 3 2]), N_g1*T, []);
        ydot_K3_g1 = reshape(ydot_g1',N_g1*T,1);
        zdot_K3_g2 = reshape(permute(zdot_g2,[1 3 2]), N_g2*T, []);
        ydot_K3_g2 = reshape(ydot_g2',N_g2*T,1);
        zdot_K3_g3 = reshape(permute(zdot_g3,[1 3 2]), N_g3*T, []);
        ydot_K3_g3 = reshape(ydot_g3',N_g3*T,1);

        % find pi         
        pi_K3_g1 = pinv(zdot_K3_g1'*zdot_K3_g1)*zdot_K3_g1'*ydot_K3_g1;
        pi_K3_g2 = pinv(zdot_K3_g2'*zdot_K3_g2)*zdot_K3_g2'*ydot_K3_g2;
        pi_K3_g3 = pinv(zdot_K3_g3'*zdot_K3_g3)*zdot_K3_g3'*ydot_K3_g3;

        % residual sum of squares
        RSS1 = (ydot_K3_g1 - zdot_K3_g1*pi_K3_g1)'*(ydot_K3_g1 - zdot_K3_g1*pi_K3_g1);
        RSS2 = (ydot_K3_g2 - zdot_K3_g2*pi_K3_g2)'*(ydot_K3_g2 - zdot_K3_g2*pi_K3_g2);
        RSS3 = (ydot_K3_g3 - zdot_K3_g3*pi_K3_g3)'*(ydot_K3_g3 - zdot_K3_g3*pi_K3_g3);

        % find sigv2hat       
        sigv2hat_K3_g1 = (1/(N_g1*(T-1)))*RSS1;
        sigv2hat_K3_g2 = (1/(N_g2*(T-1)))*RSS2;
        sigv2hat_K3_g3 = (1/(N_g3*(T-1)))*RSS3;

        % estimate standard error of sigv2hat
        vit_g1 = ydot_K3_g1 - zdot_K3_g1*pi_K3_g1;
        vit_g2 = ydot_K3_g2 - zdot_K3_g2*pi_K3_g2;
        vit_g3 = ydot_K3_g3 - zdot_K3_g3*pi_K3_g3;

        vit2bar_g1 = sum(vit_g1.^2)/(N_g1*T);
        vit2bar_g2 = sum(vit_g2.^2)/(N_g2*T);
        vit2bar_g3 = sum(vit_g3.^2)/(N_g3*T);

        varhatvit2_g1 = (1/(N_g1*(T-1)))*sum( (vit_g1.^2 - vit2bar_g1).^2 );
        varhatvit2_g2 = (1/(N_g2*(T-1)))*sum( (vit_g2.^2 - vit2bar_g2).^2 );
        varhatvit2_g3 = (1/(N_g3*(T-1)))*sum( (vit_g3.^2 - vit2bar_g3).^2 );

        SE_sigvhat_K3_g1 = ( sqrt(varhatvit2_g1)/sqrt(N_g1*T) )/(2*sqrt(sigv2hat_K3_g1));
        SE_sigvhat_K3_g2 = ( sqrt(varhatvit2_g2)/sqrt(N_g2*T) )/(2*sqrt(sigv2hat_K3_g2));
        SE_sigvhat_K3_g3 = ( sqrt(varhatvit2_g3)/sqrt(N_g3*T) )/(2*sqrt(sigv2hat_K3_g3));

        % information criterion
        % No lambda
        IC_k3_L0 = N_g1*T*log(sqrt(sigv2hat_K3_g1)) + RSS1./sigv2hat_K3_g1 ...
                 + N_g2*T*log(sqrt(sigv2hat_K3_g2)) + RSS2./sigv2hat_K3_g2 ...
                 + N_g3*T*log(sqrt(sigv2hat_K3_g3)) + RSS3./sigv2hat_K3_g3;

        IC_k3_L1 = IC_k3_L0 + lambda_s3*k;    % c_lambda=1
        
        fprintf('size of clusters if K*=3: %s\n', ...
            mat2str([size(ydot_g1,1),size(ydot_g2,1),size(ydot_g3,1)]));

    elseif k == 4 
        % construct group 1 response and regressors 
        y_g1  =  y(K_sieves==1,:);
        [N_g1,T1] = size(y_g1);
        if isnan(mLbar_S3)
            mLbar_K4g1 = floor((N_g1*T)^0.2);   % set m lower bar
        else
            mLbar_K4g1 = mLbar_S3;
        end
        ZimLbar_K4_g1 = inputs_g(y,x1,x2,x3,x4,x5,B1, B2, B3, B4, B5, B6, B7, B8, B9, B10, B11, K_sieves,1,mLbar_K4g1,Params);        
        ydot_g1 = y_g1 - mean(y_g1,2);
        zdot_g1 = ZimLbar_K4_g1 - mean(ZimLbar_K4_g1,1);
        
        % construct group 2 response and regressors
        y_g2  =  y(K_sieves==2,:);
        [N_g2,T2] = size(y_g2);
        if isnan(mLbar_S3)
            mLbar_K4g2 = floor((N_g2*T)^0.2);   % set m lower bar
        else
            mLbar_K4g2 = mLbar_S3;
        end
        ZimLbar_K4_g2 = inputs_g(y,x1,x2,x3,x4,x5,B1, B2, B3, B4, B5, B6, B7, B8, B9, B10, B11, K_sieves,2,mLbar_K4g2,Params);        
        ydot_g2 = y_g2 - mean(y_g2,2);
        zdot_g2 = ZimLbar_K4_g2 - mean(ZimLbar_K4_g2,1);
        
        % construct group 3 response and regressors
        y_g3  =  y(K_sieves==3,:);
        [N_g3,T3] = size(y_g3);
        if isnan(mLbar_S3)
            mLbar_K4g3 = floor((N_g3*T)^0.2);   % set m lower bar
        else
            mLbar_K4g3 = mLbar_S3;
        end
        ZimLbar_K4_g3 = inputs_g(y,x1,x2,x3,x4,x5,B1, B2, B3, B4, B5, B6, B7, B8, B9, B10, B11, K_sieves,3,mLbar_K4g3,Params);        
        ydot_g3 = y_g3 - mean(y_g3,2);
        zdot_g3 = ZimLbar_K4_g3 - mean(ZimLbar_K4_g3,1);
        
        % construct group 4 response and regressors
        y_g4  =  y(K_sieves==4,:);
        [N_g4,T4] = size(y_g4);
        if isnan(mLbar_S3)
            mLbar_K4g4 = floor((N_g4*T)^0.2);   % set m lower bar
        else
            mLbar_K4g4 = mLbar_S3;
        end
        ZimLbar_K4_g4 = inputs_g(y,x1,x2,x3,x4,x5,B1, B2, B3, B4, B5,B6, B7, B8, B9, B10, B11, K_sieves,4,mLbar_K4g4,Params); 

        ydot_g4 = y_g4 - mean(y_g4,2);
        zdot_g4 = ZimLbar_K4_g4 - mean(ZimLbar_K4_g4,1);
        
        % stack response and regressors into matrix form
        zdot_K4_g1 = reshape(permute(zdot_g1,[1 3 2]), N_g1*T, []);
        ydot_K4_g1 = reshape(ydot_g1',N_g1*T,1);
        zdot_K4_g2 = reshape(permute(zdot_g2,[1 3 2]), N_g2*T, []);
        ydot_K4_g2 = reshape(ydot_g2',N_g2*T,1);
        zdot_K4_g3 = reshape(permute(zdot_g3,[1 3 2]), N_g3*T, []);
        ydot_K4_g3 = reshape(ydot_g3',N_g3*T,1);
        zdot_K4_g4 = reshape(permute(zdot_g4,[1 3 2]), N_g4*T, []);
        ydot_K4_g4 = reshape(ydot_g4',N_g4*T,1);
      
        % find pihat using OLS        
        pi_K4_g1 = pinv(zdot_K4_g1'*zdot_K4_g1)*zdot_K4_g1'*ydot_K4_g1;
        pi_K4_g2 = pinv(zdot_K4_g2'*zdot_K4_g2)*zdot_K4_g2'*ydot_K4_g2;
        pi_K4_g3 = pinv(zdot_K4_g3'*zdot_K4_g3)*zdot_K4_g3'*ydot_K4_g3;
        pi_K4_g4 = pinv(zdot_K4_g4'*zdot_K4_g4)*zdot_K4_g4'*ydot_K4_g4;

        % find RSS
        RSS1 = (ydot_K4_g1 - zdot_K4_g1*pi_K4_g1)'*(ydot_K4_g1 - zdot_K4_g1*pi_K4_g1);
        RSS2 = (ydot_K4_g2 - zdot_K4_g2*pi_K4_g2)'*(ydot_K4_g2 - zdot_K4_g2*pi_K4_g2);
        RSS3 = (ydot_K4_g3 - zdot_K4_g3*pi_K4_g3)'*(ydot_K4_g3 - zdot_K4_g3*pi_K4_g3);
        RSS4 = (ydot_K4_g4 - zdot_K4_g4*pi_K4_g4)'*(ydot_K4_g4 - zdot_K4_g4*pi_K4_g4);
        
        % find sigv2hat
        sigv2hat_K4_g1 = (1/(N_g1*(T-1)))*RSS1;
        sigv2hat_K4_g2 = (1/(N_g2*(T-1)))*RSS2;
        sigv2hat_K4_g3 = (1/(N_g3*(T-1)))*RSS3;
        sigv2hat_K4_g4 = (1/(N_g4*(T-1)))*RSS4;
        
        % estimate standard error of sigv2hat
        vit_g1 = ydot_K4_g1 - zdot_K4_g1*pi_K4_g1;
        vit_g2 = ydot_K4_g2 - zdot_K4_g2*pi_K4_g2;
        vit_g3 = ydot_K4_g3 - zdot_K4_g3*pi_K4_g3;
        vit_g4 = ydot_K4_g4 - zdot_K4_g4*pi_K4_g4;

        vit2bar_g1 = sum(vit_g1.^2)/(N_g1*T);
        vit2bar_g2 = sum(vit_g2.^2)/(N_g2*T);
        vit2bar_g3 = sum(vit_g3.^2)/(N_g3*T);
        vit2bar_g4 = sum(vit_g4.^2)/(N_g4*T);

        varhatvit2_g1 = (1/(N_g1*(T-1)))*sum( (vit_g1.^2 - vit2bar_g1).^2 );
        varhatvit2_g2 = (1/(N_g2*(T-1)))*sum( (vit_g2.^2 - vit2bar_g2).^2 );
        varhatvit2_g3 = (1/(N_g3*(T-1)))*sum( (vit_g3.^2 - vit2bar_g3).^2 );
        varhatvit2_g4 = (1/(N_g4*(T-1)))*sum( (vit_g4.^2 - vit2bar_g4).^2 );

        SE_sigvhat_K4_g1 = ( sqrt(varhatvit2_g1)/sqrt(N_g1*T) )/(2*sqrt(sigv2hat_K4_g1));
        SE_sigvhat_K4_g2 = ( sqrt(varhatvit2_g2)/sqrt(N_g2*T) )/(2*sqrt(sigv2hat_K4_g2));
        SE_sigvhat_K4_g3 = ( sqrt(varhatvit2_g3)/sqrt(N_g3*T) )/(2*sqrt(sigv2hat_K4_g3));
        SE_sigvhat_K4_g4 = ( sqrt(varhatvit2_g4)/sqrt(N_g4*T) )/(2*sqrt(sigv2hat_K4_g4));

        % information criterion
        % No lambda
        IC_k4_L0 = N_g1*T*log(sqrt(sigv2hat_K4_g1)) + RSS1./sigv2hat_K4_g1 ...
                 + N_g2*T*log(sqrt(sigv2hat_K4_g2)) + RSS2./sigv2hat_K4_g2 ...
                 + N_g3*T*log(sqrt(sigv2hat_K4_g3)) + RSS3./sigv2hat_K4_g3 ...
                 + N_g4*T*log(sqrt(sigv2hat_K4_g4)) + RSS4./sigv2hat_K4_g4;

        IC_k4_L1 = IC_k4_L0 + lambda_s3*k;    % c_lambda=1
        
        fprintf('size of clusters if K*=4: %s\n', ...
            mat2str([size(ydot_g1,1),size(ydot_g2,1),size(ydot_g3,1),size(ydot_g4,1)]));

    end
end
if K_group == 1
    [C,Kstar_L1] = min(IC_k1_L1);          % single candidate -> K* = 1
    store_Kstar_L1 = Kstar_L1;
    IC0_vec = IC_k1_L0;
elseif K_group == 2
    % determine K*
    [C,Kstar_L1] = min([IC_k1_L1, IC_k2_L1]);
    store_Kstar_L1 = Kstar_L1;
    IC0_vec = [IC_k1_L0, IC_k2_L0];
elseif K_group == 3
    % determine K*
    [C,Kstar_L1] = min([IC_k1_L1, IC_k2_L1, IC_k3_L1]);
    store_Kstar_L1 = Kstar_L1;
    IC0_vec = [IC_k1_L0, IC_k2_L0, IC_k3_L0];
elseif K_group == 4
    % determine K*
    [C,Kstar_L1] = min([IC_k1_L1, IC_k2_L1, IC_k3_L1, IC_k4_L1]);
    store_Kstar_L1 = Kstar_L1;
    IC0_vec = [IC_k1_L0, IC_k2_L0, IC_k3_L0, IC_k4_L0];
end


Ksieves1 = cluster(Z_sieves,'maxclust',1);      % max number of cluster/group is ,K_group
Ksieves2 = cluster(Z_sieves,'maxclust',2);
Ksieves3 = cluster(Z_sieves,'maxclust',3);
Ksieves4 = cluster(Z_sieves,'maxclust',4);

%---------------------------------------------------------------------%
% Step 4 and 4': find c, sigu, and tau
%---------------------------------------------------------------------%    
% set option for fminsearch/fminsearchbnd
options = optimset('MaxFunEvals',1000,'MaxIter',1000,'TolFun',1.e-4,'TolX',1.e-4);

if Kstar_L1 == 1
        
    % step 4: find omegabar = [c, sigu]
    omegabar_K1 = fminsearchbnd(@(omegabar_K1) loglik_K1(omegabar_K1, y, ZimLbar_K1, ...
        sqrt(sigv2hat_K1), pi_K1), Params.omegabar0, Params.LB, Params.UB, options);
    store_omegabar_K1(1,:) = omegabar_K1;
        
    % step 4': find omegabar = [c1, sigu1, c2, sigu2, tau]
    omegabar_K1_prime = fminsearchbnd(@(omegabar_K1_prime) loglik_K1_prime(omegabar_K1_prime, y, ZimLbar_K1, ...
        sqrt(sigv2hat_K1), pi_K1), Params.omegabar0_prime, Params.LB_prime, Params.UB_prime, options);
    store_omegabar_K1_prime(1,:) = omegabar_K1_prime;
        
    % step 5': IC of latent structure for c, sigu
    for lambda_s5 = [log(N)*sqrt(N)/8 (3/4)*log(N)*sqrt(N)/8 (3/2)*log(N)*sqrt(N)/8]
        IC1_k1_prime = loglik_K1(omegabar_K1, y, ZimLbar_K1, sqrt(sigv2hat_K1), pi_K1) + lambda_s5;
        IC2_k1_prime = loglik_K1_prime(omegabar_K1_prime, y, ZimLbar_K1, sqrt(sigv2hat_K1), pi_K1) + 2*lambda_s5;

        % check if there is mixture structure in c - u_i   
        [V,MixD] = min([IC1_k1_prime,IC2_k1_prime]);

        if lambda_s5 == log(N)*sqrt(N)/8  % c_lambda = 1
            store_MixD_L1 = MixD;
        elseif lambda_s5 == (3/4)*log(N)*sqrt(N)/8 % c_lambda = 3/4
            store_MixD_L34 = MixD;
        elseif lambda_s5 == (3/2)*log(N)*sqrt(N)/8 % c_lambda = 3/2
            store_MixD_L32 = MixD;
        end
    end 
        
elseif Kstar_L1 == 2 

    % step 4: find omegabar = [c, sigu]
    omegabar_K2 = fminsearchbnd(@(omegabar_K2) loglik_K2(omegabar_K2, y, ZimLbar_K2_g1, ZimLbar_K2_g2, ...
        sqrt(sigv2hat_K2_g1), sqrt(sigv2hat_K2_g2), pi_K2_g1, pi_K2_g2, Ksieves2), Params.omegabar0, Params.LB, Params.UB, options);
    store_omegabar_K2(1,:) = omegabar_K2;

    % step 4': find omegabar = [c1, sigu1, c2, sigu2, tau]
    omegabar_K2_prime = fminsearchbnd(@(omegabar_K2_prime) loglik_K2_prime(omegabar_K2_prime, y, ZimLbar_K2_g1, ZimLbar_K2_g2, ...
        sqrt(sigv2hat_K2_g1), sqrt(sigv2hat_K2_g2), pi_K2_g1, pi_K2_g2, Ksieves2), Params.omegabar0_prime, Params.LB_prime, Params.UB_prime, options);
    store_omegabar_K2_prime(1,:) = omegabar_K2_prime;
        
    % step 5': IC of latent structure for c, sigu
    for lambda_s5 = [log(N)*sqrt(N)/8 (3/4)*log(N)*sqrt(N)/8 (3/2)*log(N)*sqrt(N)/8 ]
        IC1_k2_prime = loglik_K2(omegabar_K2, y, ZimLbar_K2_g1, ZimLbar_K2_g2, ...
            sqrt(sigv2hat_K2_g1), sqrt(sigv2hat_K2_g2), pi_K2_g1, pi_K2_g2, Ksieves2) + lambda_s5;
        IC2_k2_prime = loglik_K2_prime(omegabar_K2_prime, y, ZimLbar_K2_g1, ZimLbar_K2_g2, ...
            sqrt(sigv2hat_K2_g1), sqrt(sigv2hat_K2_g2), pi_K2_g1, pi_K2_g2, Ksieves2) + 2*lambda_s5;

        % find if there is mixture structure in c - u_i   
        [V,MixD] = min([IC1_k2_prime,IC2_k2_prime]);

        
        if lambda_s5 == log(N)*sqrt(N)/8  % c_lambda = 1
            store_MixD_L1 = MixD;
        elseif lambda_s5 == (3/4)*log(N)*sqrt(N)/8 % c_lambda = 3/4
            store_MixD_L34 = MixD;
        elseif lambda_s5 == (3/2)*log(N)*sqrt(N)/8 % c_lambda = 3/2
            store_MixD_L32 = MixD;
        end
    end 

elseif Kstar_L1 == 3

    % step 4: find omegabar = [c, sigu]
    omegabar_K3 = fminsearchbnd(@(omegabar_K3) loglik_K3(omegabar_K3, y, ZimLbar_K3_g1, ZimLbar_K3_g2, ZimLbar_K3_g3, ...
        sqrt(sigv2hat_K3_g1), sqrt(sigv2hat_K3_g2), sqrt(sigv2hat_K3_g3), ...
        pi_K3_g1, pi_K3_g2, pi_K3_g3, Ksieves3), Params.omegabar0, Params.LB, Params.UB, options);
    store_omegabar_K3(1,:) = omegabar_K3;
        
    % step 4': find omegabar = [c1, sqrt(sigu1), c2, sqrt(sigu2), tau]
    omegabar_K3_prime = fminsearchbnd(@(omegabar_K3_prime) loglik_K3_prime(omegabar_K3_prime, y, ZimLbar_K3_g1, ZimLbar_K3_g2, ZimLbar_K3_g3, ...
        sqrt(sigv2hat_K3_g1), sqrt(sigv2hat_K3_g2), sqrt(sigv2hat_K3_g3), ... 
        pi_K3_g1, pi_K3_g2, pi_K3_g3, Ksieves3), Params.omegabar0_prime, Params.LB_prime, Params.UB_prime, options);
    store_omegabar_K3_prime(1,:) = omegabar_K3_prime;
        
    % step 5': IC of latent structure for c, sigu
    for lambda_s5 = [log(N)*sqrt(N)/8 (3/4)*log(N)*sqrt(N)/8 (3/2)*log(N)*sqrt(N)/8]
        IC1_k3_prime = loglik_K3(omegabar_K3, y, ZimLbar_K3_g1, ZimLbar_K3_g2, ZimLbar_K3_g3, ...
        sqrt(sigv2hat_K3_g1), sqrt(sigv2hat_K3_g2), sqrt(sigv2hat_K3_g3), ... 
        pi_K3_g1, pi_K3_g2, pi_K3_g3, Ksieves3) + lambda_s5;

        IC2_k3_prime = loglik_K3_prime(omegabar_K3_prime, y, ZimLbar_K3_g1, ZimLbar_K3_g2, ZimLbar_K3_g3, ...
            sqrt(sigv2hat_K3_g1), sqrt(sigv2hat_K3_g2), sqrt(sigv2hat_K3_g3), ... 
            pi_K3_g1, pi_K3_g2, pi_K3_g3, Ksieves3) + 2*lambda_s5;

        % find if there is mixture structure in c - u_i  
        [V,MixD] = min([IC1_k3_prime,IC2_k3_prime]);

        
        if lambda_s5 == log(N)*sqrt(N)/8  % c_lambda = 1
            store_MixD_L1 = MixD;
        elseif lambda_s5 == (3/4)*log(N)*sqrt(N)/8 % c_lambda = 3/4
            store_MixD_L34 = MixD;
        elseif lambda_s5 == (3/2)*log(N)*sqrt(N)/8 % c_lambda = 3/2
            store_MixD_L32 = MixD;
        end
    end 

elseif Kstar_L1 == 4
        
    % step 4: find omegabar = [c, sigu]
    omegabar_K4 = fminsearchbnd(@(omegabar_K4) loglik_K4(omegabar_K4, y, ZimLbar_K4_g1, ZimLbar_K4_g2, ZimLbar_K4_g3, ZimLbar_K4_g4, ...
        sqrt(sigv2hat_K4_g1), sqrt(sigv2hat_K4_g2), sqrt(sigv2hat_K4_g3), sqrt(sigv2hat_K4_g4), ...
        pi_K4_g1, pi_K4_g2, pi_K4_g3, pi_K4_g4, Ksieves4), Params.omegabar0, Params.LB, Params.UB, options);
    store_omegabar_K4(1,:) = omegabar_K4;

    % step 4': find omegabar = [c1, sigu1, c2, sigu2, tau]
    omegabar_K4_prime = fminsearchbnd(@(omegabar_K4_prime) loglik_K4_prime(omegabar_K4_prime, y, ZimLbar_K4_g1, ZimLbar_K4_g2, ZimLbar_K4_g3, ZimLbar_K4_g4, ...
        sqrt(sigv2hat_K4_g1), sqrt(sigv2hat_K4_g2), sqrt(sigv2hat_K4_g3), sqrt(sigv2hat_K4_g4), ... 
        pi_K4_g1, pi_K4_g2, pi_K4_g3, pi_K4_g4, Ksieves4), Params.omegabar0_prime, Params.LB_prime, Params.UB_prime, options);
    store_omegabar_K4_prime(1,:) = omegabar_K4_prime;


    % step 5': IC of latent structure for c, sigu
    for lambda_s5 = [log(N)*sqrt(N)/8 (3/4)*log(N)*sqrt(N)/8 (3/2)*log(N)*sqrt(N)/8]
        IC1_k4_prime = loglik_K4(omegabar_K4, y, ZimLbar_K4_g1, ZimLbar_K4_g2, ZimLbar_K4_g3, ZimLbar_K4_g4, ...
        sqrt(sigv2hat_K4_g1), sqrt(sigv2hat_K4_g2), sqrt(sigv2hat_K4_g3), sqrt(sigv2hat_K4_g4), ... 
        pi_K4_g1, pi_K4_g2, pi_K4_g3, pi_K4_g4, Ksieves4) + lambda_s5;
        
        IC2_k4_prime = loglik_K4_prime(omegabar_K4_prime, y, ZimLbar_K4_g1, ZimLbar_K4_g2, ZimLbar_K4_g3, ZimLbar_K4_g4, ...
        sqrt(sigv2hat_K4_g1), sqrt(sigv2hat_K4_g2), sqrt(sigv2hat_K4_g3), sqrt(sigv2hat_K4_g4), ... 
        pi_K4_g1, pi_K4_g2, pi_K4_g3, pi_K4_g4, Ksieves4) + 2*lambda_s5;
            
        % find if there is mixture structure in c - u_i  
        [V,MixD] = min([IC1_k4_prime,IC2_k4_prime]);

           
        if lambda_s5 == log(N)*sqrt(N)/8  % c_lambda = 1
            store_MixD_L1 = MixD;
        elseif lambda_s5 == (3/4)*log(N)*sqrt(N)/8 % c_lambda = 3/4
            store_MixD_L34 = MixD;
        elseif lambda_s5 == (3/2)*log(N)*sqrt(N)/8 % c_lambda = 3/2
            store_MixD_L32 = MixD;
        end

    end  
end 

% optimal numer of group for frontiers
fprintf('\n optimal number of groups (K*): %d \n', Kstar_L1);
% c-u structure: =1 if unique, =2 if mixture
fprintf(' c-u structure (1=unique, 2=mixture): %d \n', store_MixD_L1);

% robustness: optimal number of groups K* under alternative STEP-3 penalty c_lambda
% K* = argmin_k ( IC_k_L0 + c_lambda*lambda_s3*k ); benchmark is c_lambda = 1.
% Printed (not saved) so the replication code stays minimal.
fprintf('\n robustness -- optimal number of groups K* under alternative step-3 penalty (c_lambda):\n');
c_lambda_grid = [3/2, 1, 3/4];
k_grid        = 1:numel(IC0_vec);
for c_lambda = c_lambda_grid
    [~, Kstar_rob] = min(IC0_vec + c_lambda*lambda_s3*k_grid);
    fprintf('   step 3 penalty term c_lambda = %5.4f :  K* = %d\n', c_lambda, Kstar_rob);
end

% robustness: c-u structure under alternative STEP-5' penalty c_tilde_lambda,
% CONDITIONAL on K* selected under the benchmark step-3 penalty (c_lambda = 1).
% c-u value: 1 = unique, 2 = mixture. Printed, not saved.
cu_lbl = {'unique', 'mixture'};
fprintf('\n robustness -- c-u structure under alternative step-5'' penalty (c_tilde_lambda),\n');
fprintf(' conditional on K* = %d selected at the benchmark step-3 penalty (c_lambda = 1):\n', Kstar_L1);
fprintf('   step 5'' penalty term c_tilde_lambda = %5.4f :  %d (%s)\n', 3/2, store_MixD_L32, cu_lbl{store_MixD_L32});
fprintf('   step 5'' penalty term c_tilde_lambda = %5.4f :  %d (%s)\n', 1,   store_MixD_L1,  cu_lbl{store_MixD_L1});
fprintf('   step 5'' penalty term c_tilde_lambda = %5.4f :  %d (%s)\n', 3/4, store_MixD_L34, cu_lbl{store_MixD_L34});

%-------------------------------------------------------------------------%
% Save frontiers
%-------------------------------------------------------------------------%
if Kstar_L1 == 1
    store_pi_K1{mLbar_K1} = pi_K1(:);
elseif Kstar_L1 == 2
    store_pi_K2_g1{mLbar_K2g1} = pi_K2_g1(:);
    store_pi_K2_g2{mLbar_K2g2} = pi_K2_g2(:);
elseif Kstar_L1 == 3
    store_pi_K3_g1{mLbar_K3g1} = pi_K3_g1(:);
    store_pi_K3_g2{mLbar_K3g2} = pi_K3_g2(:);
    store_pi_K3_g3{mLbar_K3g3} = pi_K3_g3(:);
elseif Kstar_L1 == 4
    store_pi_K4_g1{mLbar_K4g1} = pi_K4_g1(:);
    store_pi_K4_g2{mLbar_K4g2} = pi_K4_g2(:);
    store_pi_K4_g3{mLbar_K4g3} = pi_K4_g3(:);
    store_pi_K4_g4{mLbar_K4g4} = pi_K4_g4(:);
end
%-------------------------------------------------------------------------%
% Save \hat{sigma}_v
%-------------------------------------------------------------------------%
if Kstar_L1 == 1
    store_sigvhat_K1(1,:) = sqrt(sigv2hat_K1); 
elseif Kstar_L1 == 2
    store_sigvhat_K2(1,:) = sqrt([sigv2hat_K2_g1, sigv2hat_K2_g2]');
elseif Kstar_L1 == 3
    store_sigvhat_K3(1,:) = sqrt([sigv2hat_K3_g1, sigv2hat_K3_g2, sigv2hat_K3_g3]');
elseif Kstar_L1 == 4
    store_sigvhat_K4(1,:) = sqrt([sigv2hat_K4_g1, sigv2hat_K4_g2, sigv2hat_K4_g3, sigv2hat_K4_g4]');
end

fprintf('\n %%-- finished running algorithm --%%\n');

% ---------------------------------------------------------------------------
% Pack every local variable created above into the output struct.
% who is snapshotted BEFORE the helper names below exist, so they are excluded;
% eval reads each local by name. This guarantees the driver receives the full
% workspace the downstream scripts expect, without enumerating ~150 names.
% ---------------------------------------------------------------------------
local_vars = who;
out = struct();
for ii_pack = 1:numel(local_vars)
    out.(local_vars{ii_pack}) = eval(local_vars{ii_pack});
end

end
