% detector 2 pixel size = 1.4497


distSC = 1726-273;%1310;     % mm  distance between detectors
% LCLS parameters
% distance from focus to first detector = 1.4986 meters
% distance from membrane to first detector = 895 mm, between det =420 mm
% focal distance of the lenses = 317.5
% LCLS pixel size pix1 = 1.44276, pix2 = 1.44978;
% EuXFEL pixel size pi1 = 1.4338673051 pix2 = 1.20776665922
pixsize = 1.44276;             % um            % um           % pixelsize in um
pixsize2 = 1.44978;

% distance as in the manuscript of Seaberg
R1 = 478e-3;        % in meters
R2 = 1925e-3;       % in meters

maskDiam = 1500./scalef;            % diameter mask size in um for second detector

dispscale = 1;
EkeV = 9.5;
stepMot = 0.4;          % scan step of the motor in um
delta = 3.776e-6;       % Beryllium at 9.5 keV density 1.848

%-----------------------------------------------------------------------
lambda = 1.2398/EkeV.*1e-9;
kw2 = 2*pi/lambda;

% export plots and cuts 
% empty to not save anything, save folder otherwise
exportPlots =  'E:\LCLS\processing\Grating_1s'; 



%% =============================  project on the axis ======================

yper1 = m/ypeak1;       xper1 = n/xpeak1; % use period
angleRot = -atan(ypeak1/xpeak1);

dv1pd = dv1p ./(2*pi) .* ((m/(ypeak1-1))^2 + (n/(xpeak1-1))^2 ).^(1/2) .*pixsize*1e-6 ;
dv2pd = dv2p ./(2*pi) .* ((m/(ypeak1-1))^2 + (n/(xpeak1-1))^2 ).^(1/2) .*pixsize*1e-6 ;

gradH = -cos(angleRot)*dv1pd - sin(angleRot)*dv2pd;
gradV = sin(angleRot)*dv1pd - cos(angleRot)*dv2pd;
%% =============================  scale the phase ======================
% This is the formula from seaberg
gradH = gradH.*R1./R2./(R2 - R1) ./sqrt(2);  % wavefront gradient and not phase gradient
gradV = gradV.*R1./R2./(R2 - R1) ./sqrt(2);
%% Integration




figure(19)
subplot(1,3,1)
imagesc(gradV.*smmask)
ylabel('Wavefront gradient')
colorbar('SouthOutSide')
subplot(1,3,2)
imagesc(gradH.*smmask)
colorbar('SouthOutSide')
colormap gray
subplot(1,3,3)
%imagesc(z3)
colorbar('SouthOutSide')



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        integration                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mask1 = smmask;     mask2 = smmask;
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

output2 = ZernikeCalc(1:15,phi, mask1);
output2 = output2.*1e9;

exportLensFig2D


% % % %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % %                              export plots                               %
% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


return;
if ~isempty(exportPlots)
    %% export drawings
    dlmwrite(exportFileSummary,['Radius calculated with diff fitting is in meters  ' num2str((sum(Rh(smmask(:)).^Lnorme)./sum(smmask(:))).^(1/Lnorme)) '  H  ' num2str((sum(Rv(smmask(:)).^Lnorme)./sum(smmask(:))).^(1/Lnorme)) ' V' ], '-append','delimiter','');
    dlmwrite(exportFileSummary,['Standard deviation  ' num2str(std2(Rh)) '  H  ' num2str(std2(Rv)) ' V' ], '-append','delimiter','');
 
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
    
    save(fullfile(exportPlots,'phi'),'phi1','phiFull','file1','phiFull')
    
    filecuts = fullfile(exportPlots,'cuts.dat');
    dlmwrite(filecuts, 'Cut lines ','delimiter', '' );
    dlmwrite(filecuts, 'Axis (mm) - Vcut (nm) - Hcut (nm)', '-append','delimiter', '' );
    dlmwrite(filecuts, cutsout, '-append','delimiter', '\t' );
end;

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