%% this script is to compare with and without phase plate:wavefrot,
% thickness and phase

%%   load results from the two processing
noPlateG = load('E:\LCLS\processing\XSVT\noplateRescaled\environ.mat');
wiPlateG = load('E:\LCLS\processing\XSVT\plateRescaled\environ.mat');
file1 = noPlateG.file1;
absorptionnP = noPlateG.absorption;     absorptionwP = wiPlateG.absorption; 
noPlateG = noPlateG.out1;               wiPlateG = wiPlateG.out1;   
%% %%%%%%%%%%%%%%%%%%%%5 setup configuration %%%%%%%%
% detector 2 pixel size = 1.4497

dist = 420;     % mm  distance between detectors
% LCLS parameters
% distance from focus to first detector = 1.4986 meters
% distance from membrane to first detector = 895 mm
% distance between the detectors = 420
% focal distance of the lenses = 317.5
% LCLS pixel size pix1 = 1.44276, pix2 = ;or 1.6105
% EuXFEL pixel size pi1 = 1.4338673051 pix2 = 1.20776665922
interpsecim = 2.62/2;%(1.4338673051/1.20776665922)^(1);  
MagnLensToDet = 317/1.4986; % We use this factor to rescale the wavefront from the first detector position to the exit of the lens

wolverine = 99;                          % did Wolverine played with your data? filter a keep only the low frequencies of the gradients
pixsize = 1.44276;                      % um            % um           % pixelsize in um
pixsize2 = 1.44978.*interpsecim;        % because you interpolated in the first script
maskDiam = 1150;                        % diameter mask size in um

EkeV = 9.5;
delta = 3.776e-6;%
lambda = 1.2398/EkeV .*1e-9;
kw2 = 2*pi/lambda; 
% export plots and cuts 
% empty to not save anything, save folder otherwise
exportPlots = '/data/bm05/imaging/seb/LCLS/processing/XSVT/comparison'; 
saveEnvironment = '/data/bm05/imaging/seb/LCLS/processing/XSVT/comparison'; 
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           build masks                                   $
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('Center','var'),     Center = [];        end
[smmask, bigmask,Center] = maskbuilder(absorptionnP,maskDiam/pixsize,Center);
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        correction                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

out1 = noPlateG - wiPlateG;

% load('/data/bm05/imaging/seb/xfel/trackvectOneLensXFELlarge steps_2.mat')
% load('/data/bm05/imaging/seb/xfel/trackvectNoLensXFELlarge_binning.mat');
% erosion to remove shot noise - shouldnt be more than a couuple of %
vert  = myerosion(out1(:,:,1).*bigmask,3,3);
horiz = myerosion(out1(:,:,2).*bigmask,3,3);

[JH, JV] = meshgrid(1:size(vert,2),1:(size(vert,1)));
%put pixel size difference into 
horiz = horiz  - JH.*(pixsize2 - pixsize)./pixsize2;    
vert  = vert   - JV.*(pixsize2 - pixsize)./pixsize2;

psb = cat(3,horiz,vert);
vert  = vert(min(JV(bigmask==1)):max(JV(bigmask==1)) , min(JH(bigmask==1)):max(JH(bigmask==1)),:);
horiz = horiz(min(JV(bigmask==1)):max(JV(bigmask==1)),min(JH(bigmask==1)):max(JH(bigmask==1)),:);

% remove detector rotation
[vp,vr] = legxy(vert,smmask);
[hp,hr] = legxy(horiz,smmask);

% and put in rad by the same line
vert2   = (vert -  vr) .* pixsize2./dist.*1e-3;
horiz2  = (horiz - hr) .* pixsize2./dist.*1e-3;
% remove tilt of the beam
horiz2 = horiz2 - mean(horiz2(:));
vert2  = vert2  - mean(vert2(:));
ps = cat(3,horiz2,vert2);
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        integration                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mask1 = smmask;     mask2 = bigmask;

% pixsize = pixsize.* MagnLensToDet;
% ps = ps./MagnLensToDet;
%calculate gradient error
[p1,p2] = gradient_error(ps(:,:,1),ps(:,:,1),mask1);
[p3,p4] = gradient_error(ps(:,:,2),ps(:,:,2),mask1);
    
% calculate differential focal length
raddiff1 = 1/(mean(mean(diff(p1,1,2))./(pixsize*1e-6)));
raddiff2 = 1/(mean(mean(diff(p4,1,1))./(pixsize*1e-6)));
disp(['Differential focal length with diff fitting is in meters  ' num2str(raddiff1) '  H  ' num2str(raddiff2) ' V ']);
disp(['For a mask of diameter = ' num2str(round(maskDiam)) 'um']);

% calculate a rmse moyenne
dd = ((mean(mean(diff(p1,1,2)))./2 + mean(mean(diff(p4,1,1)))./2 ));
p11 = p1./mean(mean(diff(p1,1,2))).*dd;
p44 = p4./mean(mean(diff(p4,1,1))).*dd;
raddiff11 = 1/(mean(mean(diff(p11,1,2)./(pixsize*1e-6))));
raddiff22 = 1/(mean(mean(diff(p44,1,1)./(pixsize*1e-6))));

%extract wavefront gradient error and put in urad (with astigmatism)

xg = (ps(:,:,1)-p11).*1e6;    xg = xg - mean(xg(mask1));
yg = (ps(:,:,2)-p44).*1e6;    yg = yg - mean(yg(mask1));
if wolverine ~= 0
    % wolverine effect removal
    [output1,~] = ZernikeCalc(1:wolverine,xg, mask1,'STANDARD');
    xg = sum(output1,3);
    [output1,~] = ZernikeCalc(1:wolverine,yg, mask1,'STANDARD');
    yg = sum(output1,3);
end 

% do integration to get wavefront in meter
phi = WftSolveLSChol(yg.*mask1,xg.*mask1,pixsize,mask1).*1e-12;
phi = phi- mean(phi(mask1(:)));% remove piston

%extract wavefront gradient error and put in urad (without astigmatism and defocus)
xg1 = (ps(:,:,1)-p2 -p1).*1e6;    xg1 = xg1 - mean(xg1(mask1));
yg1 = (ps(:,:,2)-p3 -p4).*1e6;    yg1 = yg1 - mean(yg1(mask1));
if wolverine ~= 0
    % wolverine effect removal
    [output1,~] = ZernikeCalc(1:wolverine,xg1, mask1,'STANDARD');
    xg1 = sum(output1,3);
    [output1,~] = ZernikeCalc(1:wolverine,yg1, mask1,'STANDARD');
    yg1 = sum(output1,3);
end

phi1 = WftSolveLSChol(yg1.*mask1,xg1.*mask1,pixsize,mask1).*1e-12;
phi1 = phi1- mean(phi1(mask1(:)));% remove pistonfull map phase
phiFull = WftSolveLSChol(( vert2 - mean( vert2(mask1))).*mask1,(horiz2 - mean(horiz2(mask1))).*mask1,pixsize,ones(size(horiz2))).*1e-12;
phiFull = phiFull- mean(phiFull(mask1(:)));% remove piston

[p1,p2] = gradient_error(phi1,phi1,mask1);
phi1 = phi1 - p1 - p2;

%thickness depends on delta value
thickErr = phi./delta.*1e6;% put in microns by the way
thickErr = thickErr - mean(thickErr(mask1(:)));% remove mean
%thickness error without astimatism
thickErr1 = phi1./delta.*1e6;% put in microns by the way
thickErr1 = thickErr1 - mean(thickErr1(mask1(:)));% remove mean


phiD = phi1;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               Without plate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

out1 = noPlateG ;

% load('/data/bm05/imaging/seb/xfel/trackvectOneLensXFELlarge steps_2.mat')
% load('/data/bm05/imaging/seb/xfel/trackvectNoLensXFELlarge_binning.mat');
% erosion to remove shot noise - shouldnt be more than a couuple of %
vert  = myerosion(out1(:,:,1).*bigmask,3,3);
horiz = myerosion(out1(:,:,2).*bigmask,3,3);

[JH, JV] = meshgrid(1:size(vert,2),1:(size(vert,1)));
%put pixel size difference into 
horiz = horiz  - JH.*(pixsize2 - pixsize)./pixsize2;    
vert  = vert   - JV.*(pixsize2 - pixsize)./pixsize2;

psb = cat(3,horiz,vert);
vert  = vert(min(JV(bigmask==1)):max(JV(bigmask==1)) , min(JH(bigmask==1)):max(JH(bigmask==1)),:);
horiz = horiz(min(JV(bigmask==1)):max(JV(bigmask==1)),min(JH(bigmask==1)):max(JH(bigmask==1)),:);

% remove detector rotation
[vp,vr] = legxy(vert,smmask);
[hp,hr] = legxy(horiz,smmask);

% and put in rad by the same line
vert2   = (vert -  vr) .* pixsize2./dist.*1e-3;
horiz2  = (horiz - hr) .* pixsize2./dist.*1e-3;
% remove tilt of the beam
horiz2 = horiz2 - mean(horiz2(:));
vert2  = vert2  - mean(vert2(:));
ps = cat(3,horiz2,vert2);
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        integration                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mask1 = smmask;     mask2 = bigmask;

% pixsize = pixsize.* MagnLensToDet;
% ps = ps./MagnLensToDet;
%calculate gradient error
[p1,p2] = gradient_error(ps(:,:,1),ps(:,:,1),mask1);
[p3,p4] = gradient_error(ps(:,:,2),ps(:,:,2),mask1);
    
% calculate differential focal length
raddiff1 = 1/(mean(mean(diff(p1,1,2))./(pixsize*1e-6)));
raddiff2 = 1/(mean(mean(diff(p4,1,1))./(pixsize*1e-6)));
disp(['Differential focal length with diff fitting is in meters  ' num2str(raddiff1) '  H  ' num2str(raddiff2) ' V ']);
disp(['For a mask of diameter = ' num2str(round(maskDiam)) 'um']);

% calculate a rmse moyenne
dd = ((mean(mean(diff(p1,1,2)))./2 + mean(mean(diff(p4,1,1)))./2 ));
p11 = p1./mean(mean(diff(p1,1,2))).*dd;
p44 = p4./mean(mean(diff(p4,1,1))).*dd;
raddiff11 = 1/(mean(mean(diff(p11,1,2)./(pixsize*1e-6))));
raddiff22 = 1/(mean(mean(diff(p44,1,1)./(pixsize*1e-6))));

%extract wavefront gradient error and put in urad (with astigmatism)

xg = (ps(:,:,1)-p11).*1e6;    xg = xg - mean(xg(mask1));
yg = (ps(:,:,2)-p44).*1e6;    yg = yg - mean(yg(mask1));
if wolverine ~= 0
    % wolverine effect removal
    [output1,~] = ZernikeCalc(1:wolverine,xg, mask1,'STANDARD');
    xg = sum(output1,3);
    [output1,~] = ZernikeCalc(1:wolverine,yg, mask1,'STANDARD');
    yg = sum(output1,3);
end 

% do integration to get wavefront in meter
phi = WftSolveLSChol(yg.*mask1,xg.*mask1,pixsize,mask1).*1e-12;
phi = phi- mean(phi(mask1(:)));% remove piston

%extract wavefront gradient error and put in urad (without astigmatism and defocus)
xg1 = (ps(:,:,1)-p2 -p1).*1e6;    xg1 = xg1 - mean(xg1(mask1));
yg1 = (ps(:,:,2)-p3 -p4).*1e6;    yg1 = yg1 - mean(yg1(mask1));
if wolverine ~= 0
    % wolverine effect removal
    [output1,~] = ZernikeCalc(1:wolverine,xg1, mask1,'STANDARD');
    xg1 = sum(output1,3);
    [output1,~] = ZernikeCalc(1:wolverine,yg1, mask1,'STANDARD');
    yg1 = sum(output1,3);
end

phi1 = WftSolveLSChol(yg1.*mask1,xg1.*mask1,pixsize,mask1).*1e-12;
phi1 = phi1- mean(phi1(mask1(:)));% remove pistonfull map phase
phiFull = WftSolveLSChol(( vert2 - mean( vert2(mask1))).*mask1,(horiz2 - mean(horiz2(mask1))).*mask1,pixsize,ones(size(horiz2))).*1e-12;
phiFull = phiFull- mean(phiFull(mask1(:)));% remove piston

[p1,p2] = gradient_error(phi1,phi1,mask1);
phi1 = phi1 - p1 - p2;

%thickness depends on delta value
thickErr = phi./delta.*1e6;% put in microns by the way
thickErr = thickErr - mean(thickErr(mask1(:)));% remove mean
%thickness error without astimatism
thickErr1 = phi1./delta.*1e6;% put in microns by the way
thickErr1 = thickErr1 - mean(thickErr1(mask1(:)));% remove mean


phiWo = phi1;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               With plate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
out1 = wiPlateG;

% load('/data/bm05/imaging/seb/xfel/trackvectOneLensXFELlarge steps_2.mat')
% load('/data/bm05/imaging/seb/xfel/trackvectNoLensXFELlarge_binning.mat');
% erosion to remove shot noise - shouldnt be more than a couuple of %
vert  = myerosion(out1(:,:,1).*bigmask,3,3);
horiz = myerosion(out1(:,:,2).*bigmask,3,3);

[JH, JV] = meshgrid(1:size(vert,2),1:(size(vert,1)));
%put pixel size difference into 
horiz = horiz  - JH.*(pixsize2 - pixsize)./pixsize2;    
vert  = vert   - JV.*(pixsize2 - pixsize)./pixsize2;

psb = cat(3,horiz,vert);
vert  = vert(min(JV(bigmask==1)):max(JV(bigmask==1)) , min(JH(bigmask==1)):max(JH(bigmask==1)),:);
horiz = horiz(min(JV(bigmask==1)):max(JV(bigmask==1)),min(JH(bigmask==1)):max(JH(bigmask==1)),:);

% remove detector rotation
[vp,vr] = legxy(vert,smmask);
[hp,hr] = legxy(horiz,smmask);

% and put in rad by the same line
vert2   = (vert -  vr) .* pixsize2./dist.*1e-3;
horiz2  = (horiz - hr) .* pixsize2./dist.*1e-3;
% remove tilt of the beam
horiz2 = horiz2 - mean(horiz2(:));
vert2  = vert2  - mean(vert2(:));
ps = cat(3,horiz2,vert2);
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        integration                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mask1 = smmask;     mask2 = bigmask;

% pixsize = pixsize.* MagnLensToDet;
% ps = ps./MagnLensToDet;
%calculate gradient error
[p1,p2] = gradient_error(ps(:,:,1),ps(:,:,1),mask1);
[p3,p4] = gradient_error(ps(:,:,2),ps(:,:,2),mask1);
    
% calculate differential focal length
raddiff1 = 1/(mean(mean(diff(p1,1,2))./(pixsize*1e-6)));
raddiff2 = 1/(mean(mean(diff(p4,1,1))./(pixsize*1e-6)));
disp(['Differential focal length with diff fitting is in meters  ' num2str(raddiff1) '  H  ' num2str(raddiff2) ' V ']);
disp(['For a mask of diameter = ' num2str(round(maskDiam)) 'um']);

% calculate a rmse moyenne
dd = ((mean(mean(diff(p1,1,2)))./2 + mean(mean(diff(p4,1,1)))./2 ));
p11 = p1./mean(mean(diff(p1,1,2))).*dd;
p44 = p4./mean(mean(diff(p4,1,1))).*dd;
raddiff11 = 1/(mean(mean(diff(p11,1,2)./(pixsize*1e-6))));
raddiff22 = 1/(mean(mean(diff(p44,1,1)./(pixsize*1e-6))));

%extract wavefront gradient error and put in urad (with astigmatism)

xg = (ps(:,:,1)-p11).*1e6;    xg = xg - mean(xg(mask1));
yg = (ps(:,:,2)-p44).*1e6;    yg = yg - mean(yg(mask1));
if wolverine ~= 0
    % wolverine effect removal
    [output1,~] = ZernikeCalc(1:wolverine,xg, mask1,'STANDARD');
    xg = sum(output1,3);
    [output1,~] = ZernikeCalc(1:wolverine,yg, mask1,'STANDARD');
    yg = sum(output1,3);
end 

% do integration to get wavefront in meter
phi = WftSolveLSChol(yg.*mask1,xg.*mask1,pixsize,mask1).*1e-12;
phi = phi- mean(phi(mask1(:)));% remove piston

%extract wavefront gradient error and put in urad (without astigmatism and defocus)
xg1 = (ps(:,:,1)-p2 -p1).*1e6;    xg1 = xg1 - mean(xg1(mask1));
yg1 = (ps(:,:,2)-p3 -p4).*1e6;    yg1 = yg1 - mean(yg1(mask1));
if wolverine ~= 0
    % wolverine effect removal
    [output1,~] = ZernikeCalc(1:wolverine,xg1, mask1,'STANDARD');
    xg1 = sum(output1,3);
    [output1,~] = ZernikeCalc(1:wolverine,yg1, mask1,'STANDARD');
    yg1 = sum(output1,3);
end

phi1 = WftSolveLSChol(yg1.*mask1,xg1.*mask1,pixsize,mask1).*1e-12;
phi1 = phi1- mean(phi1(mask1(:)));% remove pistonfull map phase
phiFull = WftSolveLSChol(( vert2 - mean( vert2(mask1))).*mask1,(horiz2 - mean(horiz2(mask1))).*mask1,pixsize,ones(size(horiz2))).*1e-12;
phiFull = phiFull- mean(phiFull(mask1(:)));% remove piston

[p1,p2] = gradient_error(phi1,phi1,mask1);
phi1 = phi1 - p1 - p2;

%thickness depends on delta value
thickErr = phi./delta.*1e6;% put in microns by the way
thickErr = thickErr - mean(thickErr(mask1(:)));% remove mean
%thickness error without astimatism
thickErr1 = phi1./delta.*1e6;% put in microns by the way
thickErr1 = thickErr1 - mean(thickErr1(mask1(:)));% remove mean


phiW = phi1;


exportFig2D_wWo