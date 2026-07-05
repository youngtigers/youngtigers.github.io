function [B0, B1, B2, B3, B4, B5, B6, B7, B8, B9, B10, B11] = GenBasis(Params)
T  = Params.T;
t  = (1:T)';
B0  = ones(T,1);
B1  = sqrt(2)*cos(    pi*t/T);
B2  = sqrt(2)*cos(  2*pi*t/T);
B3  = sqrt(2)*cos(  3*pi*t/T);
B4  = sqrt(2)*cos(  4*pi*t/T);
B5  = sqrt(2)*cos(  5*pi*t/T);
B6  = sqrt(2)*cos(  6*pi*t/T);
B7  = sqrt(2)*cos(  7*pi*t/T);
B8  = sqrt(2)*cos(  8*pi*t/T);
B9  = sqrt(2)*cos(  9*pi*t/T);
B10 = sqrt(2)*cos( 10*pi*t/T);
B11 = sqrt(2)*cos( 11*pi*t/T);
end
