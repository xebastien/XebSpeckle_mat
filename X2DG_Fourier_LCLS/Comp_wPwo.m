% ===============================================================
pixreach.V  = 3; % number of interpixel distance to correlate
pixreach.H  = 3;
scanDimOrd  = 1; % )How was made the mesh scan: 1st dim is vert or horizontal ? boolean 1 or 0
PitchG      = 4;     StepS = 0.4;    nStep = 26;
padFact     = 5;
scalef      = 1;
maskDiam    = 1550./scalef; %1650            % diameter mask size in um

pixsize = 1.44978;
delta = 3.776e-6;       % Beryllium at 9.5 keV density 1.848
EkeV = 9.5;
lambda = 1.2398/EkeV.*1e-9;
kw2 = 2*pi/lambda;
% distance as in the manuscript of Seaberg
R1 = 478e-3;        % in meters
R2 = 1925e-3;       % in meters

magnif = 11;

per_diag = sqrt(2)*magnif*PitchG;
PeriodScan = StepS*nStep/PitchG;

saveName = 'E:\LCLS\processing\Grating_scanning\Checkboard_Plate.mat';

%% ======================================================================
% Compare with and without plate wavefront of the beam from metrology using
% the scanning Fourier garting processing

load('E:\LCLS\processing\Grating_scanning\Checkboard_Plate.mat')
dvwp1 = dv1;     dvwp2 = dv2;   wpmask = smmask;


load('E:\LCLS\processing\Grating_scanning\Checkboard_noPlate.mat')
dvw1 = dv1;     dvw2 = -dv2;
angleRot = -atan(pi/4);
%% =============================  project/rotate on the axis ======================

[reconst1, reconst2] = gradient_error(dvw1,dvw1,smmask);
dvp1 = dvw1 - reconst1 - reconst2;
dvp1 = dvp1 - mean(dvp1(smmask(:)));
[reconst1, reconst2] = gradient_error(dvw2,dvw2,smmask);
dvp2 = dvw2 - reconst1 - reconst2;
dv2p = dvp2 - mean(dvp2(smmask(:)));


dvp1 = dvp1 ./(2*pi) .*    PitchG .*sqrt(2) *1e-6 ; % scale with the period of the peak in Fourier
dvp2 = dvp2 ./(2*pi) .*    PitchG .*sqrt(2) *1e-6;

gradH = (-cos(angleRot)*dvp1 - sin(angleRot) *dvp2);  % rotation
gradV = (sin(angleRot) *dvp1 - cos(angleRot) *dvp2);
%% =============================  scale the phase ======================
% This is the formula from seaberg
gradH = gradH .*R1 ./R2 .*(R2 - R1);  % wavefront gradient and not phase gradient
gradV = gradV .*R1 ./R2 .*(R2 - R1);


%% %%%%%%%%%%%%%%%%%%%%%%%  Integrating with plate %%%%%%%%%%%%%%%%%%%%

mask1 = wpmask;     mask2 = wpmask;
%extract wavefront gradient error and put in urad (with astigmatism)
vert2 = gradV;      horiz2 = gradH;
raddiff1 = [];      raddiff2 = [];

xg = gradH.*1e6;    xg = xg - mean(xg(mask1));
yg = gradV.*1e6;    yg = yg - mean(yg(mask1));
% do integration to get wavefront in meter
phi = WftSolveLSChol(yg.*mask1,xg.*mask1,pixsize,mask1).*1e-12;
phi1 = phi- mean(phi(mask1(:)));% remove piston
phi  = phi1;    phiFull = phi;

%thickness depends on delta value
thickErr  = phi./delta.*1e6;% put in microns by the way
thickErr  = thickErr - mean(thickErr(mask1(:)));% remove mean
%thickness error without astimatism
thickErr1 = phi1./delta.*1e6;% put in microns by the way
thickErr1 = thickErr1 - mean(thickErr1(mask1(:)));% remove mean


phiWo = phi1;

%% =============================  project/rotate on the axis ======================

[reconst1, reconst2] = gradient_error(dvwp1,dvwp1,smmask);
dvp1 = dvwp1 - reconst1 - reconst2;
dvp1 = dvp1 - mean(dvp1(smmask(:)));
[reconst1, reconst2] = gradient_error(dvwp2,dvwp2,smmask);
dvp2 = dvwp2 - reconst1 - reconst2;
dv2p = dvp2 - mean(dvp2(smmask(:)));


dvp1 = dvp1 ./(2*pi) .*    PitchG .*sqrt(2) *1e-6 ; % scale with the period of the peak in Fourier
dvp2 = dvp2 ./(2*pi) .*    PitchG .*sqrt(2) *1e-6;

gradH = (-cos(angleRot)*dvp1 - sin(angleRot) *dvp2);  % rotation
gradV = (sin(angleRot) *dvp1 - cos(angleRot) *dvp2);
%% =============================  scale the phase ======================
% This is the formula from seaberg
gradH = gradH .*R1 ./R2 .*(R2 - R1);  % wavefront gradient and not phase gradient
gradV = gradV .*R1 ./R2 .*(R2 - R1);


%% %%%%%%%%%%%%%%%%%%%%%%%  Integrating with plate %%%%%%%%%%%%%%%%%%%%

mask1 = wpmask;     mask2 = wpmask;
%extract wavefront gradient error and put in urad (with astigmatism)
vert2 = gradV; horiz2 = gradH;
raddiff1 = [];raddiff2 = [];

xg = gradH.*1e6;    xg = xg - mean(xg(mask1));
yg = gradV.*1e6;    yg = yg - mean(yg(mask1));
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


phiD = phiWo - phiW ;


exportFig2D_wWo