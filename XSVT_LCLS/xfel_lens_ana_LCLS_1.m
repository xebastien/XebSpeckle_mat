% detector 2 pixel size = 1.4497

dist = 420;     % mm  distance between detectors
% LCLS parameters
% distance from focus to first detector = 1.4986 meters
% distance from membrane to first detector = 895 mm
% distance between the detectors = 420
% focal distance of the lenses = 317.5
% LCLS pixel size pix1 = 1.44276, pix2 = ;or 1.6105
% EuXFEL pixel size pi1 = 1.4338673051 pix2 = 1.20776665922

MagnLensToDet = 317/1.4986; % We use this factor to rescale the wavefront from the first detector position to the exit of the lens

wolverine = 1;                          % did Wolverine played with your data? filter a keep only the low frequencies of the gradients
pixsize = 1.44276;                      % um            % um           % pixelsize in um
pixsize2 = 1.44978.*interpsecim;        % because you interpolated in the first script
maskDiam = 1200;                        % diameter mask size in um

EkeV = 9.5;
delta = 3.776e-6;%
lambda = 1.2398/EkeV .*1e-9;
kw2 = 2*pi/lambda; 

%videoName = '~/LCLS_XSTpulses.avi';
% export plots and cuts 
% empty to not save anything, save folder otherwise
exportPlots = '/data/bm05/imaging/seb/LCLS/processing/XSVT/noplateRescaled'; 
saveEnvironment = '/data/bm05/imaging/seb/LCLS/processing/XSVT/noplateRescaled'; 
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           for the records                               $
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isempty(exportPlots)
    exportFileSummary = fullfile(exportPlots,'parameters_processing.txt');

    dlmwrite(exportFileSummary, 'Processing parameters for second step','delimiter', '' );
    dlmwrite(exportFileSummary, ['pixsize2 = ' num2str(pixsize2) ' um'], '-append','delimiter', '' );
    dlmwrite(exportFileSummary, ['pixsize1 = ' num2str(pixsize) ' um'], '-append','delimiter', '' );
    dlmwrite(exportFileSummary, ['dist = ' num2str(dist) ' mm'], '-append','delimiter', '' );
    dlmwrite(exportFileSummary, ['maskDiam = ' num2str(maskDiam) ' um'], '-append','delimiter', '' );
    dlmwrite(exportFileSummary, ['Energy = ' num2str(EkeV) ' keV'], '-append','delimiter', '' );
    dlmwrite(exportFileSummary, ['delta = ' num2str(delta)], '-append','delimiter', '' );
    dlmwrite(exportFileSummary, ['lambda = ' num2str(lambda) ' m'], '-append','delimiter','');
    dlmwrite(exportFileSummary, ['winsize = ' num2str(winsize)], '-append','delimiter','');
    dlmwrite(exportFileSummary, ['winvect = ' num2str(winvect)], '-append','delimiter','');
    dlmwrite(exportFileSummary, ['ROI1 = ' num2str(ROI1)], '-append','delimiter','');
    %dlmwrite(exportFileSummary, 'Processed:', '-append','delimiter', '' );
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           build masks                                   $
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('Center','var'),     Center = [];        end
[smmask, bigmask,Center] = maskbuilder(absorption,maskDiam/pixsize,Center);
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
if wolverine
    % wolverine effect removal
    [output1,~] = ZernikeCalc(1:99,xg, mask1,'STANDARD');
    xg = sum(output1,3);
    [output1,~] = ZernikeCalc(1:99,yg, mask1,'STANDARD');
    yg = sum(output1,3);
end 

% do integration to get wavefront in meter
phi = WftSolveLSChol(yg.*mask1,xg.*mask1,pixsize,mask1).*1e-12;
phi = phi- mean(phi(mask1(:)));% remove piston

%extract wavefront gradient error and put in urad (without astigmatism and defocus)
xg1 = (ps(:,:,1)-p2 -p1).*1e6;    xg1 = xg1 - mean(xg1(mask1));
yg1 = (ps(:,:,2)-p3 -p4).*1e6;    yg1 = yg1 - mean(yg1(mask1));
if wolverine
    % wolverine effect removal
    [output1,~] = ZernikeCalc(1:99,xg1, mask1,'STANDARD');
    xg1 = sum(output1,3);
    [output1,~] = ZernikeCalc(1:99,yg1, mask1,'STANDARD');
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


exportLensFig2D


% % % %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % %                              export plots                               %
% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



if ~isempty(exportPlots)
    %% export drawings
    dlmwrite(exportFileSummary,['Radius calculated with diff fitting is in meters  ' num2str(raddiff1) '  H  ' num2str(raddiff2) ' V '], '-append','delimiter','');
    %dlmwrite(exportFileSummary,['Standard deviation  ' num2str(std2(Rh)) '  H  ' num2str(std2(Rv)) ' V' ], '-append','delimiter','');
 
    print(1,'-dpdf',fullfile(exportPlots,'Wavefront_gradient.pdf'));
    saveas(1, fullfile(exportPlots,'Wavefront_gradient.svg'));
    print(1,'-dpng',fullfile(exportPlots,'Wavefront_gradient.png'));
    
    print(2,'-dpdf',fullfile(exportPlots,'Wavefront_gradient_error.pdf'));
    saveas(2, fullfile(exportPlots,'Wavefront_gradient_error.svg'));
    print(2,'-dpng',fullfile(exportPlots,'Wavefront_gradient_error.png'));

    print(3,'-dpdf',fullfile(exportPlots,'Maks.pdf'));
    saveas(3, fullfile(exportPlots,'Maks.svg'));
    print(3,'-dpng',fullfile(exportPlots,'Maks.png'));
    
    print(4,'-dpdf',fullfile(exportPlots,'Wavefront_error.pdf'));
    saveas(4, fullfile(exportPlots,'Wavefront_error.svg'));
    print(4,'-dpng',fullfile(exportPlots,'Wavefront_error.png'));
    
    %print(5,'-dpdf',fullfile(exportPlots,'Thickness_3D_error.pdf'));
    saveas(5, fullfile(exportPlots,'Wavefront_gradient_error.svg'));
    print(5,'-dpng',fullfile(exportPlots,'Thickness_3D_error.png'));
    
    print(6,'-dpdf',fullfile(exportPlots,'WavefrontError_cuts.pdf'));
    saveas(6, fullfile(exportPlots,'WavefrontError_cuts.svg'));
    print(6,'-dpng',fullfile(exportPlots,'WavefrontError_cuts.png'));
    
    save(fullfile(exportPlots,'phi'),'phi1','phiFull','file1','phiFull','out1')
    
    filecuts = fullfile(exportPlots,'cuts.dat');
    dlmwrite(filecuts, 'Cut lines ','delimiter', '' );
    dlmwrite(filecuts, 'Axis (mm) - Vcut (nm) - Hcut (nm)', '-append','delimiter', '' );
    dlmwrite(filecuts, cutsout, '-append','delimiter', '\t' );
end
if ~isempty(saveEnvironment)
    save(fullfile(saveEnvironment,'environ.mat'),'out1','absorption','file1');
end
    
    

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