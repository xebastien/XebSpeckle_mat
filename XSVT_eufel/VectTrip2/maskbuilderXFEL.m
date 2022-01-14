function [smmask, bigmask,Center] = maskbuilderXFEL(absorption,maskSize,CenterIN)


[m,n] = size(absorption);
% BUILD A MASK
[X1,Y1] = meshgrid(1:n,1:m);
X = X1 - round(n/2);  Y = Y1-round(m/2);
mask0 = (X.^2+Y.^2).^(1/2);
exponent = exp(-(((X).^2 + (Y).^2)./(2*(maskSize./2)^2)));

exponent(mask0 < (maskSize./2)) = 1;
mask3 = (mask0 < (maskSize./2.*1.3)).*exponent;


a = sort(reshape(absorption,1,[]),'descend');
maskthreshold = sum(mask3(:) ==1)./numel(mask3);

b = a(round(length(a)*maskthreshold));
d1 = (absorption > b) ;

[x_offset,y_offset] = fshift(d1,mask3);    
Center(2) = x_offset + n/2;   Center(1) = y_offset + m/2;
if nargin ==3,
    Center = CenterIN;
end;

%  mask
[X1,Y1] = meshgrid(1:n,1:m);

X = X1 - Center(2);  Y = Y1-Center(1);
mask0 = (X.^2+Y.^2).^(1/2);

mask2 = (mask0 < (maskSize./2));
bigmask = mask2 == 1;


smmask = bigmask(min(Y1(bigmask==1)):max(Y1(bigmask==1)),min(X1(bigmask==1)):max(X1(bigmask==1)));


