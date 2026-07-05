function Q = inputs_g(y,x1,x2,x3,x4,x5,B1,B2,B3,B4,B5,B6,B7,B8,B9,B10,B11,K_sieves,g,mLbar,Params)
x1_g = x1(K_sieves==g,:);
x2_g = x2(K_sieves==g,:);
x3_g = x3(K_sieves==g,:);
x4_g = x4(K_sieves==g,:);
x5_g = x5(K_sieves==g,:);
y_g  =  y(K_sieves==g,:);
[N_g, T]  = size(y_g);
B0        = ones(T,1);
BmLbar    = [B0, B1, B2, B3, B4, B5, B6, B7, B8, B9, B10, B11];
zBm_g     = NaN(T, mLbar*5, N_g);
ZimLbar_g = NaN(T, mLbar*6-1, N_g);
for i = 1:N_g
    for t = 1:T
        b = BmLbar(t, 1:mLbar);
        zBm_g(t,:,i)     = [x1_g(i,t)*b, x2_g(i,t)*b, x3_g(i,t)*b, x4_g(i,t)*b, x5_g(i,t)*b];
        ZimLbar_g(t,:,i) = [BmLbar(t,2:mLbar), zBm_g(t,:,i)];
    end
end
Q = ZimLbar_g;
end
