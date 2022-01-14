%% Compare with and without phase plate beam phase
% Built figure and export
%
% For processing grating single shot  

load('E:\LCLS\processing\Grating_1s\Checkboard_Plate.mat')
wpGradH = -gradH;   wpGradV = -gradV;  wpmask = smmask; wpCenter = Center;


load('E:\LCLS\processing\Grating_1s\Checkboard_noPlate.mat')
npGradH = -gradH;   npGradV = -gradV; npmask = smmask;  npCenter = Center;


pixsize = 1.44978;
delta = 3.776e-6;       % Beryllium at 9.5 keV density 1.848
EkeV = 9.5;
lambda = 1.2398/EkeV.*1e-9;
kw2 = 2*pi/lambda;

%% %%%%%%%%%%%%%%%%%%%%%%%  Integrating with plate %%%%%%%%%%%%%%%%%%%%

mask1 = wpmask;     mask2 = wpmask;
%extract wavefront gradient error and put in urad (with astigmatism)

xg = wpGradH.*1e6;    xg = xg - mean(xg(mask1));
yg = wpGradV.*1e6;    yg = yg - mean(yg(mask1));
% do integration to get wavefront in meter
phi = WftSolveLSChol(yg.*mask1,xg.*mask1,pixsize,mask1).*1e-12;
phi1 = phi- mean(phi(mask1(:)));% remove piston
phi = phi1;phiFull = phi;

%thickness depends on delta value
thickErr = phi./delta.*1e6;% put in microns by the way
thickErr = thickErr - mean(thickErr(mask1(:)));% remove mean
%thickness error without astimatism
thickErr1 = phi1./delta.*1e6;% put in microns by the way
thickErr1 = thickErr1 - mean(thickErr1(mask1(:)));% remove mean


phiW = phi1;


%% %%%%%%%%%%%%%%%%%%%%%%%  Integrating with plate %%%%%%%%%%%%%%%%%%%%

mask1 = npmask;     mask2 = npmask;
%extract wavefront gradient error and put in urad (with astigmatism)

xg = npGradH.*1e6;    xg = xg - mean(xg(mask1));
yg = npGradV.*1e6;    yg = yg - mean(yg(mask1));
% do integration to get wavefront in meter
phi = WftSolveLSChol(yg.*mask1,xg.*mask1,pixsize,mask1).*1e-12;
phi1 = phi- mean(phi(mask1(:)));% remove piston
phi = phi1;phiFull = phi;

%thickness depends on delta value
thickErr = phi./delta.*1e6;% put in microns by the way
thickErr = thickErr - mean(thickErr(mask1(:)));% remove mean
%thickness error without astimatism
thickErr1 = phi1./delta.*1e6;% put in microns by the way
thickErr1 = thickErr1 - mean(thickErr1(mask1(:)));% remove mean


phiWo = phi1;



%% %%%%%%%%%%%%%%%%%%%%%%%  Integrating with plate %%%%%%%%%%%%%%%%%%%%

mask1 = npmask;     mask2 = npmask;
%extract wavefront gradient error and put in urad (with astigmatism)
vert2 = npGradV - wpGradV; horiz2 = npGradH - wpGradH;
raddiff1 = [];raddiff2 = [];

xg = (npGradH - wpGradH).*1e6;    xg = xg - mean(xg(mask1));
yg = (npGradV - wpGradV).*1e6;    yg = yg - mean(yg(mask1));
% do integration to get wavefront in meter
phi = WftSolveLSChol(yg.*mask1,xg.*mask1,pixsize,mask1).*1e-12;
phi1 = phi- mean(phi(mask1(:)));% remove piston
phi = phi1;phiFull = phi;

%thickness depends on delta value
thickErr = phi./delta.*1e6;% put in microns by the way
thickErr = thickErr - mean(thickErr(mask1(:)));% remove mean
%thickness error without astimatism
thickErr1 = phi1./delta.*1e6;% put in microns by the way
thickErr1 = thickErr1 - mean(thickErr1(mask1(:)));% remove mean


phiD = phi1;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%    Plots   %%%%%%%%%%%%%%%%%%%%%%%%%%%
phiWo = phiWo - mean(phiWo(:));
phiW = phiW - mean(phiW(:));
phiD = phiD - mean(phiD(:));

exportFig2D_wWo