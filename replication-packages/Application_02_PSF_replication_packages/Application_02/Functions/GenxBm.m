function xBm = GenxBm(x1,x2,x3,x4,x5,B1,B2,B3,B4,B5,B6,B7,B8,B9,B10,B11,Params)
N    = Params.N;
T    = Params.T;
m    = Params.m;
B0   = ones(T,1);
Ball = [B0, B1, B2, B3, B4, B5, B6, B7, B8, B9, B10, B11];
xBm  = NaN(T, m*5, N);
for i = 1:N
    for t = 1:T
        b = Ball(t, 1:m);
        xBm(t,:,i) = [x1(i,t)*b, x2(i,t)*b, x3(i,t)*b, x4(i,t)*b, x5(i,t)*b];
    end
end
end
