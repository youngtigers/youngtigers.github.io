%% SE for \sigma_{v}
% Computed inline in main_application.m at the point of estimation (step 1),
% immediately after each sigv2hat estimate:
%   K* = 1: SE_sigvhat_K1            (main_application.m ~line 192)
%   K* = 2: SE_sigvhat_K2_g1/g2      (main_application.m ~lines 260-261)
%   K* = 3: SE_sigvhat_K3_g1..g3     (main_application.m ~lines 346-348)
%   K* = 4: SE_sigvhat_K4_g1..g4     (main_application.m ~lines 455-458)
% Formula: SE(sigvhat) = ( sqrt(Var(v_it^2))/sqrt(N*T) ) / (2*sqrt(sigv2hat)),
% i.e. SE of the variance sigv2hat times the delta-method factor 1/(2*sigvhat).
% Nothing is computed in this section; see those lines.


%% SE for frontiers
for t = 1:T
    BmLbar(t,:) = [B0(t), B1(t), B2(t), B3(t), B4(t), B5(t), B6(t), B7(t), B8(t), B9(t), B10(t), B11(t)]'; % T x m lower bar
end

if Kstar_L1 == 1

% storage
MB_K1 = NaN(p+1,mLbar_K1-1+mLbar_K1*p,T);
S_K1 = NaN(p+1,p+1,T);


Q_K1_zz = zdotK1'*zdotK1/(N*T);

for t = 1:T
    MB_K1(:,:,t) = [BmLbar(t,2:mLbar_K1), zeros(1,mLbar_K1),    zeros(1,mLbar_K1),    zeros(1,mLbar_K1),    zeros(1,mLbar_K1),    zeros(1,mLbar_K1);
                    zeros(1,mLbar_K1-1),  BmLbar(t,1:mLbar_K1), zeros(1,mLbar_K1),    zeros(1,mLbar_K1),    zeros(1,mLbar_K1),    zeros(1,mLbar_K1);
                    zeros(1,mLbar_K1-1),  zeros(1,mLbar_K1),    BmLbar(t,1:mLbar_K1), zeros(1,mLbar_K1),    zeros(1,mLbar_K1),    zeros(1,mLbar_K1);
                    zeros(1,mLbar_K1-1),  zeros(1,mLbar_K1),    zeros(1,mLbar_K1),    BmLbar(t,1:mLbar_K1), zeros(1,mLbar_K1),    zeros(1,mLbar_K1);
                    zeros(1,mLbar_K1-1),  zeros(1,mLbar_K1),    zeros(1,mLbar_K1),    zeros(1,mLbar_K1),    BmLbar(t,1:mLbar_K1), zeros(1,mLbar_K1);
                    zeros(1,mLbar_K1-1),  zeros(1,mLbar_K1),    zeros(1,mLbar_K1),    zeros(1,mLbar_K1),    zeros(1,mLbar_K1),    BmLbar(t,1:mLbar_K1)];
    S_K1(:,:,t) = sigv2hat_K1*MB_K1(:,:,t)*pinv(Q_K1_zz)*MB_K1(:,:,t)'/mLbar_K1;
end


elseif Kstar_L1 == 2
N_K2_g1 = size(ZimLbar_K2_g1,3);
N_K2_g2 = size(ZimLbar_K2_g2,3);

% storage
MB_K2_g1 = NaN(p+1,mLbar_K2g1-1+mLbar_K2g1*p,T);
MB_K2_g2 = NaN(p+1,mLbar_K2g2-1+mLbar_K2g2*p,T);
S_K2_g1 = NaN(p+1,p+1,T);
S_K2_g2 = NaN(p+1,p+1,T);

Q_K2_g1_zz = zdot_K2_g1'*zdot_K2_g1/(N_K2_g1*T);
Q_K2_g2_zz = zdot_K2_g2'*zdot_K2_g2/(N_K2_g2*T);

for t = 1:T
    MB_K2_g1(:,:,t) = [BmLbar(t,2:mLbar_K2g1), zeros(1,mLbar_K2g1),    zeros(1,mLbar_K2g1),    zeros(1,mLbar_K2g1),    zeros(1,mLbar_K2g1),    zeros(1,mLbar_K2g1);
                       zeros(1,mLbar_K2g1-1),  BmLbar(t,1:mLbar_K2g1), zeros(1,mLbar_K2g1),    zeros(1,mLbar_K2g1),    zeros(1,mLbar_K2g1),    zeros(1,mLbar_K2g1);
                       zeros(1,mLbar_K2g1-1),  zeros(1,mLbar_K2g1),    BmLbar(t,1:mLbar_K2g1), zeros(1,mLbar_K2g1),    zeros(1,mLbar_K2g1),    zeros(1,mLbar_K2g1);
                       zeros(1,mLbar_K2g1-1),  zeros(1,mLbar_K2g1),    zeros(1,mLbar_K2g1),    BmLbar(t,1:mLbar_K2g1), zeros(1,mLbar_K2g1),    zeros(1,mLbar_K2g1);
                       zeros(1,mLbar_K2g1-1),  zeros(1,mLbar_K2g1),    zeros(1,mLbar_K2g1),    zeros(1,mLbar_K2g1),    BmLbar(t,1:mLbar_K2g1), zeros(1,mLbar_K2g1);
                       zeros(1,mLbar_K2g1-1),  zeros(1,mLbar_K2g1),    zeros(1,mLbar_K2g1),    zeros(1,mLbar_K2g1),    zeros(1,mLbar_K2g1),    BmLbar(t,1:mLbar_K2g1)];
    S_K2_g1(:,:,t) = sigv2hat_K2_g1*MB_K2_g1(:,:,t)*pinv(Q_K2_g1_zz)*MB_K2_g1(:,:,t)'/mLbar_K2g1;

    MB_K2_g2(:,:,t) = [BmLbar(t,2:mLbar_K2g2), zeros(1,mLbar_K2g2),    zeros(1,mLbar_K2g2),    zeros(1,mLbar_K2g2),    zeros(1,mLbar_K2g2),    zeros(1,mLbar_K2g2);
                       zeros(1,mLbar_K2g2-1),  BmLbar(t,1:mLbar_K2g2), zeros(1,mLbar_K2g2),    zeros(1,mLbar_K2g2),    zeros(1,mLbar_K2g2),    zeros(1,mLbar_K2g2);
                       zeros(1,mLbar_K2g2-1),  zeros(1,mLbar_K2g2),    BmLbar(t,1:mLbar_K2g2), zeros(1,mLbar_K2g2),    zeros(1,mLbar_K2g2),    zeros(1,mLbar_K2g2);
                       zeros(1,mLbar_K2g2-1),  zeros(1,mLbar_K2g2),    zeros(1,mLbar_K2g2),    BmLbar(t,1:mLbar_K2g2), zeros(1,mLbar_K2g2),    zeros(1,mLbar_K2g2);
                       zeros(1,mLbar_K2g2-1),  zeros(1,mLbar_K2g2),    zeros(1,mLbar_K2g2),    zeros(1,mLbar_K2g2),    BmLbar(t,1:mLbar_K2g2), zeros(1,mLbar_K2g2);
                       zeros(1,mLbar_K2g2-1),  zeros(1,mLbar_K2g2),    zeros(1,mLbar_K2g2),    zeros(1,mLbar_K2g2),    zeros(1,mLbar_K2g2),    BmLbar(t,1:mLbar_K2g2)];

    S_K2_g2(:,:,t) = sigv2hat_K2_g2*MB_K2_g2(:,:,t)*pinv(Q_K2_g2_zz)*MB_K2_g2(:,:,t)'/mLbar_K2g2;
end

elseif Kstar_L1 == 3
N_K3_g1 = size(ZimLbar_K3_g1,3);
N_K3_g2 = size(ZimLbar_K3_g2,3);
N_K3_g3 = size(ZimLbar_K3_g3,3);

Q_K3_g1_zz = zdot_K3_g1'*zdot_K3_g1/(N_K3_g1*T);
Q_K3_g2_zz = zdot_K3_g2'*zdot_K3_g2/(N_K3_g2*T);
Q_K3_g3_zz = zdot_K3_g3'*zdot_K3_g3/(N_K3_g3*T);

% storage
MB_K3_g1 = NaN(p+1,mLbar_K3g1-1+mLbar_K3g1*p,T);
MB_K3_g2 = NaN(p+1,mLbar_K3g2-1+mLbar_K3g2*p,T);
MB_K3_g3 = NaN(p+1,mLbar_K3g3-1+mLbar_K3g3*p,T);
S_K3_g1 = NaN(p+1,p+1,T);
S_K3_g2 = NaN(p+1,p+1,T);
S_K3_g3 = NaN(p+1,p+1,T);

for t = 1:T
    MB_K3_g1(:,:,t) = [BmLbar(t,2:mLbar_K3g1), zeros(1,mLbar_K3g1),    zeros(1,mLbar_K3g1),    zeros(1,mLbar_K3g1),    zeros(1,mLbar_K3g1),    zeros(1,mLbar_K3g1);
                       zeros(1,mLbar_K3g1-1),  BmLbar(t,1:mLbar_K3g1), zeros(1,mLbar_K3g1),    zeros(1,mLbar_K3g1),    zeros(1,mLbar_K3g1),    zeros(1,mLbar_K3g1);
                       zeros(1,mLbar_K3g1-1),  zeros(1,mLbar_K3g1),    BmLbar(t,1:mLbar_K3g1), zeros(1,mLbar_K3g1),    zeros(1,mLbar_K3g1),    zeros(1,mLbar_K3g1);
                       zeros(1,mLbar_K3g1-1),  zeros(1,mLbar_K3g1),    zeros(1,mLbar_K3g1),    BmLbar(t,1:mLbar_K3g1), zeros(1,mLbar_K3g1),    zeros(1,mLbar_K3g1);
                       zeros(1,mLbar_K3g1-1),  zeros(1,mLbar_K3g1),    zeros(1,mLbar_K3g1),    zeros(1,mLbar_K3g1),    BmLbar(t,1:mLbar_K3g1), zeros(1,mLbar_K3g1);
                       zeros(1,mLbar_K3g1-1),  zeros(1,mLbar_K3g1),    zeros(1,mLbar_K3g1),    zeros(1,mLbar_K3g1),    zeros(1,mLbar_K3g1),    BmLbar(t,1:mLbar_K3g1)];

    S_K3_g1(:,:,t) = sigv2hat_K3_g1*MB_K3_g1(:,:,t)*pinv(Q_K3_g1_zz)*MB_K3_g1(:,:,t)'/mLbar_K3g1;

    MB_K3_g2(:,:,t) = [BmLbar(t,2:mLbar_K3g2), zeros(1,mLbar_K3g2),    zeros(1,mLbar_K3g2),    zeros(1,mLbar_K3g2),    zeros(1,mLbar_K3g2),    zeros(1,mLbar_K3g2);
                       zeros(1,mLbar_K3g2-1),  BmLbar(t,1:mLbar_K3g2), zeros(1,mLbar_K3g2),    zeros(1,mLbar_K3g2),    zeros(1,mLbar_K3g2),    zeros(1,mLbar_K3g2);
                       zeros(1,mLbar_K3g2-1),  zeros(1,mLbar_K3g2),    BmLbar(t,1:mLbar_K3g2), zeros(1,mLbar_K3g2),    zeros(1,mLbar_K3g2),    zeros(1,mLbar_K3g2);
                       zeros(1,mLbar_K3g2-1),  zeros(1,mLbar_K3g2),    zeros(1,mLbar_K3g2),    BmLbar(t,1:mLbar_K3g2), zeros(1,mLbar_K3g2),    zeros(1,mLbar_K3g2);
                       zeros(1,mLbar_K3g2-1),  zeros(1,mLbar_K3g2),    zeros(1,mLbar_K3g2),    zeros(1,mLbar_K3g2),    BmLbar(t,1:mLbar_K3g2), zeros(1,mLbar_K3g2);
                       zeros(1,mLbar_K3g2-1),  zeros(1,mLbar_K3g2),    zeros(1,mLbar_K3g2),    zeros(1,mLbar_K3g2),    zeros(1,mLbar_K3g2),    BmLbar(t,1:mLbar_K3g2)];

    S_K3_g2(:,:,t) = sigv2hat_K3_g2*MB_K3_g2(:,:,t)*pinv(Q_K3_g2_zz)*MB_K3_g2(:,:,t)'/mLbar_K3g2;

    MB_K3_g3(:,:,t) = [BmLbar(t,2:mLbar_K3g3), zeros(1,mLbar_K3g3),    zeros(1,mLbar_K3g3),    zeros(1,mLbar_K3g3),    zeros(1,mLbar_K3g3),    zeros(1,mLbar_K3g3);
                       zeros(1,mLbar_K3g3-1),  BmLbar(t,1:mLbar_K3g3), zeros(1,mLbar_K3g3),    zeros(1,mLbar_K3g3),    zeros(1,mLbar_K3g3),    zeros(1,mLbar_K3g3);
                       zeros(1,mLbar_K3g3-1),  zeros(1,mLbar_K3g3),    BmLbar(t,1:mLbar_K3g3), zeros(1,mLbar_K3g3),    zeros(1,mLbar_K3g3),    zeros(1,mLbar_K3g3);
                       zeros(1,mLbar_K3g3-1),  zeros(1,mLbar_K3g3),    zeros(1,mLbar_K3g3),    BmLbar(t,1:mLbar_K3g3), zeros(1,mLbar_K3g3),    zeros(1,mLbar_K3g3);
                       zeros(1,mLbar_K3g3-1),  zeros(1,mLbar_K3g3),    zeros(1,mLbar_K3g3),    zeros(1,mLbar_K3g3),    BmLbar(t,1:mLbar_K3g3), zeros(1,mLbar_K3g3);
                       zeros(1,mLbar_K3g3-1),  zeros(1,mLbar_K3g3),    zeros(1,mLbar_K3g3),    zeros(1,mLbar_K3g3),    zeros(1,mLbar_K3g3),    BmLbar(t,1:mLbar_K3g3)];

    S_K3_g3(:,:,t) = sigv2hat_K3_g3*MB_K3_g3(:,:,t)*pinv(Q_K3_g3_zz)*MB_K3_g3(:,:,t)'/mLbar_K3g3;
end

elseif Kstar_L1 == 4
N_K4_g1 = size(ZimLbar_K4_g1,3);
N_K4_g2 = size(ZimLbar_K4_g2,3);
N_K4_g3 = size(ZimLbar_K4_g3,3);
N_K4_g4 = size(ZimLbar_K4_g4,3);

Q_K4_g1_zz = zdot_K4_g1'*zdot_K4_g1/(N_K4_g1*T);
Q_K4_g2_zz = zdot_K4_g2'*zdot_K4_g2/(N_K4_g2*T);
Q_K4_g3_zz = zdot_K4_g3'*zdot_K4_g3/(N_K4_g3*T);
Q_K4_g4_zz = zdot_K4_g4'*zdot_K4_g4/(N_K4_g4*T);

% storage
MB_K4_g1 = NaN(p+1,mLbar_K4g1-1+mLbar_K4g1*p,T);
MB_K4_g2 = NaN(p+1,mLbar_K4g2-1+mLbar_K4g2*p,T);
MB_K4_g3 = NaN(p+1,mLbar_K4g3-1+mLbar_K4g3*p,T);
MB_K4_g4 = NaN(p+1,mLbar_K4g4-1+mLbar_K4g4*p,T);
S_K4_g1 = NaN(p+1,p+1,T);
S_K4_g2 = NaN(p+1,p+1,T);
S_K4_g3 = NaN(p+1,p+1,T);
S_K4_g4 = NaN(p+1,p+1,T);

for t = 1:T
    MB_K4_g1(:,:,t) = [BmLbar(t,2:mLbar_K4g1), zeros(1,mLbar_K4g1),    zeros(1,mLbar_K4g1),    zeros(1,mLbar_K4g1),    zeros(1,mLbar_K4g1),    zeros(1,mLbar_K4g1);
                       zeros(1,mLbar_K4g1-1),  BmLbar(t,1:mLbar_K4g1), zeros(1,mLbar_K4g1),    zeros(1,mLbar_K4g1),    zeros(1,mLbar_K4g1),    zeros(1,mLbar_K4g1);
                       zeros(1,mLbar_K4g1-1),  zeros(1,mLbar_K4g1),    BmLbar(t,1:mLbar_K4g1), zeros(1,mLbar_K4g1),    zeros(1,mLbar_K4g1),    zeros(1,mLbar_K4g1);
                       zeros(1,mLbar_K4g1-1),  zeros(1,mLbar_K4g1),    zeros(1,mLbar_K4g1),    BmLbar(t,1:mLbar_K4g1), zeros(1,mLbar_K4g1),    zeros(1,mLbar_K4g1);
                       zeros(1,mLbar_K4g1-1),  zeros(1,mLbar_K4g1),    zeros(1,mLbar_K4g1),    zeros(1,mLbar_K4g1),    BmLbar(t,1:mLbar_K4g1), zeros(1,mLbar_K4g1);
                       zeros(1,mLbar_K4g1-1),  zeros(1,mLbar_K4g1),    zeros(1,mLbar_K4g1),    zeros(1,mLbar_K4g1),    zeros(1,mLbar_K4g1),    BmLbar(t,1:mLbar_K4g1)];

    S_K4_g1(:,:,t) = sigv2hat_K4_g1*MB_K4_g1(:,:,t)*pinv(Q_K4_g1_zz)*MB_K4_g1(:,:,t)'/mLbar_K4g1;

    MB_K4_g2(:,:,t) = [BmLbar(t,2:mLbar_K4g2), zeros(1,mLbar_K4g2),    zeros(1,mLbar_K4g2),    zeros(1,mLbar_K4g2),    zeros(1,mLbar_K4g2),    zeros(1,mLbar_K4g2);
                       zeros(1,mLbar_K4g2-1),  BmLbar(t,1:mLbar_K4g2), zeros(1,mLbar_K4g2),    zeros(1,mLbar_K4g2),    zeros(1,mLbar_K4g2),    zeros(1,mLbar_K4g2);
                       zeros(1,mLbar_K4g2-1),  zeros(1,mLbar_K4g2),    BmLbar(t,1:mLbar_K4g2), zeros(1,mLbar_K4g2),    zeros(1,mLbar_K4g2),    zeros(1,mLbar_K4g2);
                       zeros(1,mLbar_K4g2-1),  zeros(1,mLbar_K4g2),    zeros(1,mLbar_K4g2),    BmLbar(t,1:mLbar_K4g2), zeros(1,mLbar_K4g2),    zeros(1,mLbar_K4g2);
                       zeros(1,mLbar_K4g2-1),  zeros(1,mLbar_K4g2),    zeros(1,mLbar_K4g2),    zeros(1,mLbar_K4g2),    BmLbar(t,1:mLbar_K4g2), zeros(1,mLbar_K4g2);
                       zeros(1,mLbar_K4g2-1),  zeros(1,mLbar_K4g2),    zeros(1,mLbar_K4g2),    zeros(1,mLbar_K4g2),    zeros(1,mLbar_K4g2),    BmLbar(t,1:mLbar_K4g2)];

    S_K4_g2(:,:,t) = sigv2hat_K4_g2*MB_K4_g2(:,:,t)*pinv(Q_K4_g2_zz)*MB_K4_g2(:,:,t)'/mLbar_K4g2;

    MB_K4_g3(:,:,t) = [BmLbar(t,2:mLbar_K4g3), zeros(1,mLbar_K4g3),    zeros(1,mLbar_K4g3),    zeros(1,mLbar_K4g3),    zeros(1,mLbar_K4g3),    zeros(1,mLbar_K4g3);
                       zeros(1,mLbar_K4g3-1),  BmLbar(t,1:mLbar_K4g3), zeros(1,mLbar_K4g3),    zeros(1,mLbar_K4g3),    zeros(1,mLbar_K4g3),    zeros(1,mLbar_K4g3);
                       zeros(1,mLbar_K4g3-1),  zeros(1,mLbar_K4g3),    BmLbar(t,1:mLbar_K4g3), zeros(1,mLbar_K4g3),    zeros(1,mLbar_K4g3),    zeros(1,mLbar_K4g3);
                       zeros(1,mLbar_K4g3-1),  zeros(1,mLbar_K4g3),    zeros(1,mLbar_K4g3),    BmLbar(t,1:mLbar_K4g3), zeros(1,mLbar_K4g3),    zeros(1,mLbar_K4g3);
                       zeros(1,mLbar_K4g3-1),  zeros(1,mLbar_K4g3),    zeros(1,mLbar_K4g3),    zeros(1,mLbar_K4g3),    BmLbar(t,1:mLbar_K4g3), zeros(1,mLbar_K4g3);
                       zeros(1,mLbar_K4g3-1),  zeros(1,mLbar_K4g3),    zeros(1,mLbar_K4g3),    zeros(1,mLbar_K4g3),    zeros(1,mLbar_K4g3),    BmLbar(t,1:mLbar_K4g3)];

    S_K4_g3(:,:,t) = sigv2hat_K4_g3*MB_K4_g3(:,:,t)*pinv(Q_K4_g3_zz)*MB_K4_g3(:,:,t)'/mLbar_K4g3;

    MB_K4_g4(:,:,t) = [BmLbar(t,2:mLbar_K4g4), zeros(1,mLbar_K4g4),    zeros(1,mLbar_K4g4),    zeros(1,mLbar_K4g4),    zeros(1,mLbar_K4g4),    zeros(1,mLbar_K4g4);
                       zeros(1,mLbar_K4g4-1),  BmLbar(t,1:mLbar_K4g4), zeros(1,mLbar_K4g4),    zeros(1,mLbar_K4g4),    zeros(1,mLbar_K4g4),    zeros(1,mLbar_K4g4);
                       zeros(1,mLbar_K4g4-1),  zeros(1,mLbar_K4g4),    BmLbar(t,1:mLbar_K4g4), zeros(1,mLbar_K4g4),    zeros(1,mLbar_K4g4),    zeros(1,mLbar_K4g4);
                       zeros(1,mLbar_K4g4-1),  zeros(1,mLbar_K4g4),    zeros(1,mLbar_K4g4),    BmLbar(t,1:mLbar_K4g4), zeros(1,mLbar_K4g4),    zeros(1,mLbar_K4g4);
                       zeros(1,mLbar_K4g4-1),  zeros(1,mLbar_K4g4),    zeros(1,mLbar_K4g4),    zeros(1,mLbar_K4g4),    BmLbar(t,1:mLbar_K4g4), zeros(1,mLbar_K4g4);
                       zeros(1,mLbar_K4g4-1),  zeros(1,mLbar_K4g4),    zeros(1,mLbar_K4g4),    zeros(1,mLbar_K4g4),    zeros(1,mLbar_K4g4),    BmLbar(t,1:mLbar_K4g4)];

    S_K4_g4(:,:,t) = sigv2hat_K4_g4*MB_K4_g4(:,:,t)*pinv(Q_K4_g4_zz)*MB_K4_g4(:,:,t)'/mLbar_K4g4;
end
end
