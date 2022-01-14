%% Load and compare wavefront with and without phase plate


pixsize = 1.44978;
delta = 3.776e-6;       % Beryllium at 9.5 keV density 1.848
EkeV = 9.5;
lambda = 1.2398/EkeV .*1e-9;
kw2 = 2*pi/lambda;




load('E:\LCLS\processing\XSS2D\plateSpeckle2D_det2_mat')
phiW = -phi1./2;
Center=[];
maskDiam = 1500;


load('E:\LCLS\processing\XSS2D\noplateSpeckle2D_det2_mat')
if ~exist('Center','var'),     Center = [];        end;
[smmask, bigmask,Center] = maskbuilder(file1,maskDiam/pixsize+18,Center);

phiWo = -phi1.*3;

phiD = phiWo - phiW; phi = phiD;

mask1 = smmask;

exportFig2D_wWo