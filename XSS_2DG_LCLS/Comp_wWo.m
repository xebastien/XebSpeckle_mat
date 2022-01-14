%% Load and compare wavefront with and without phase plate


pixsize = 1.44978;
delta = 3.776e-6;       % Beryllium at 9.5 keV density 1.848
EkeV = 9.5;
lambda = 1.2398/EkeV .*1e-9;
kw2 = 2*pi/lambda;


load('E:\LCLS\processing\XSS2DG\plateDiamondmat.mat')
phiW = -phi1.*2;


load('E:\LCLS\processing\XSS2DG\noplateDiamondmat.mat')
phiWo = -phi1.*2;

phiD = phiWo - phiW; phi = phiD;

mask1 = smmask;

exportFig2D_wWo