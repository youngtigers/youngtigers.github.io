function Q = inputs_g0(y,x1,x2,x3,x4,x5,B1,B2,B3,B4,B5,B6,B7,B8,B9,B10,B11,mLbar,Params)
[N,T]   = size(y);
B0      = ones(T,1);
BmLbar  = [B0, B1, B2, B3, B4, B5, B6, B7, B8, B9, B10, B11];
zBm     = NaN(T, mLbar*5, N);
ZimLbar = NaN(T, mLbar*6-1, N);
for i = 1:N
    for t = 1:T
        b = BmLbar(t, 1:mLbar);
        zBm(t,:,i)     = [x1(i,t)*b, x2(i,t)*b, x3(i,t)*b, x4(i,t)*b, x5(i,t)*b];
        ZimLbar(t,:,i) = [BmLbar(t,2:mLbar), zBm(t,:,i)];
    end
end
Q = ZimLbar;
end
