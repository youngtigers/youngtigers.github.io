function Bm = GenBm(B1,B2,B3,B4,B5,B6,B7,B8,B9,B10,B11,Params)
m  = Params.m;
B  = [B1, B2, B3, B4, B5, B6, B7, B8, B9, B10, B11];
Bm = B(:, 1:m-1);
end
