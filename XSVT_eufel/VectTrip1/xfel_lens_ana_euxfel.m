% detector 2 pixel size = 1.4497

load('/data/bm05/imaging/seb/EU-XFEL/processing/XSVT_oneLens/output1.mat','file1','out1','bigmask')



interpsecim = 1.25;%(1.4338673051/1.2077666592
dist = 470;     % mm  distance between detectors
% LCLS parameters
% distance from focus to first detector = 1.4986 meters
% distance from membrane to first detector = 895 mm
% focal distance of the lenses = 317.5
pixsize = 1.20776665922;%       1.44276;             % um            % um           % pixelsize in um
pixsize2 = 1.44978./interpsecim;
maskDiam = 870;            % diameter mask size in um
Energy = 8.27;
delta = 3.776e-6;%
lambda = 1.2398/Energy.*1e-9;
kw2 = 2*pi/lambda;



CenterIN = [556 565];
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           build masks                                   $
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[smmask, bigmask,Center] = maskbuilderXFEL(file1,maskDiam/pixsize,CenterIN);
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        correction                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%calculate gradient error
[p1,p2] = gradient_error(ps(:,:,1),ps(:,:,1),mask1);
[p3,p4] = gradient_error(ps(:,:,2),ps(:,:,2),mask1);
    
% calculate differential focal length
raddiff1 = 1/(mean2(diff(p1,1,2)./(pixsize*1e-6)));
raddiff2 = 1/(mean2(diff(p4,1,1)./(pixsize*1e-6)));
disp(['Differential focal length with diff fitting is in meters  ' num2str(raddiff1) '  H  ' num2str(raddiff2) ' V ']);
disp(['For a mask of diameter = ' num2str(round(maskDiam)) 'um']);

% calculate a rmse moyenne
dd = ((mean2(diff(p1,1,2))./2 + mean2(diff(p4,1,1))./2 ));
p11 = p1./mean2(diff(p1,1,2)).*dd;
p44 = p4./mean2(diff(p4,1,1)).*dd;
raddiff11 = 1/(mean2(diff(p11,1,2)./(pixsize*1e-6)));
raddiff22 = 1/(mean2(diff(p44,1,1)./(pixsize*1e-6)));

%extract wavefront gradient error and put in urad (with astigmatism)

xg = (ps(:,:,1)-p11).*1e6;    xg = xg - mean(xg(mask1));
yg = (ps(:,:,2)-p44).*1e6;    yg = yg - mean(yg(mask1));
% do integration to get wavefront in meter
phi = WftSolveLSChol(yg.*mask1,xg.*mask1,pixsize,mask1).*1e-12;
phi = phi- mean(phi(mask1(:)));% remove piston

%extract wavefront gradient error and put in urad (without astigmatism and defocus)
xg1 = (ps(:,:,1)-p2 -p1).*1e6;    xg1 = xg1 - mean(xg1(mask1));
yg1 = (ps(:,:,2)-p3 -p4).*1e6;    yg1 = yg1 - mean(yg1(mask1));
phi1 = WftSolveLSChol(yg1.*mask1,xg1.*mask1,pixsize,mask1).*1e-12;
phi1 = phi1- mean(phi1(mask1(:)));% remove pistonfull map phase
phiFull = WftSolveLSChol(( vert2 - mean( vert2(mask1))).*mask1,(horiz2 - mean(horiz2(mask1))).*mask1,pixsize,ones(size(horiz2))).*1e-12;
phiFull = phiFull- mean(phiFull(mask1(:)));% remove piston

%thickness depends on delta value
thickErr = phi./delta.*1e6;% put in microns by the way
thickErr = thickErr - mean(thickErr(mask1(:)));% remove mean
%thickness error without astimatism
thickErr1 = phi1./delta.*1e6;% put in microns by the way
thickErr1 = thickErr1 - mean(thickErr1(mask1(:)));% remove mean

output2 = ZernikeCalc(1:15,phi, mask1);
output2 = output2.*1e9;

exportLensFig2D
return;
% % % %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % %                               plots                                     %
% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % barcolor = [1 0 0].*0.63;
% % 
% % [m1,n1] = size(horiz2);
% % m = (1:1:m1).*pixsize./1000;
% % n = (1:1:n1).*pixsize./1000;
% % 
% % figure(30)
% % subplot(2,2,1)
% % imagesc(n,m,vert2.*1e6 ,[-1 1].*500)
% % set(gca,'DataAspectRatio',[1 1 1]')
% % colormap gray
% % %colormap(flipud(CM1))
% % set(gca,'YTick',0);set(gca,'XTick',0);set(gca,'TickLength',[0 0]);
% % title('Vertical wavefront gradient')
% % t = colorbar('EastOutside','peer',gca);
% % set(get(t,'ylabel'),'String', 'urad','FontName','TimesNewRoman','FontSize',14);
% % scalebar
% % scalebar('Location','southeast','Colour',barcolor,'Bold','True','Unit','mm','ScaleLength',round(n(end))./4)
% % 
% % subplot(2,2,2)
% % imagesc(n,m,horiz2.*1e6 ,[-1 1].*500)
% % set(gca,'DataAspectRatio',[1 1 1]')
% % colormap gray
% % %colormap(flipud(CM1))
% % set(gca,'YTick',0);set(gca,'XTick',0);set(gca,'TickLength',[0 0]);
% % title('Horizontal wavefront gradient')
% % t = colorbar('EastOutside','peer',gca);
% % set(get(t,'ylabel'),'String', 'urad','FontName','TimesNewRoman','FontSize',14);
% % scalebar
% % scalebar('Location','southeast','Colour',barcolor,'Bold','True','Unit','mm','ScaleLength',round(n(end))./4)
% % 
% % subplot(2,2,3)
% % imagesc(n,m,phi5.*1e9.*bigmask)
% % set(gca,'DataAspectRatio',[1 1 1]')
% % colormap gray
% % %colormap(flipud(CM1))
% % set(gca,'YTick',0);set(gca,'XTick',0);set(gca,'TickLength',[0 0]);
% % title('Wavefront')
% % t = colorbar('EastOutside','peer',gca);
% % set(get(t,'ylabel'),'String', 'nm','FontName','TimesNewRoman','FontSize',14);
% % scalebar
% % scalebar('Location','southeast','Colour',barcolor,'Bold','True','Unit','mm','ScaleLength',round(n(end))./4)
% % % xlabel('mm');ylabel('mm')
% % % colorbar
% % phi5 = phiFull - min(phiFull(mask1));
% % phi5(~mask1) = nan;
% % subplot(2,2,4)
% % surf(phi5.*1e9)
% % xlabel('mm');ylabel('mm');zlabel('nm');
% % shading interp
% % title('Wavefront')