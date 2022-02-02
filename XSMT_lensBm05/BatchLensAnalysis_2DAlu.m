% batch processing for lens.
% calculate the wavefront gradient from speckle vector tracking
% for 2D lenses
% astigmatism considered
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                          
% ------------------folders paths----------------%
rootFolder = '/data/bm05/inhouse/thomas/180308lenses/processing/';
% ASlenses (1D)   Elenses (2D) Flenses (2D)  MNlenses (1D)  WLlenses  (1D)
% ------------------export folder----------------%
% check that you have the writing permission in this folder
exportFolder = rootFolder;%'/data/bm05/inhouse/thomas/170410lenstest1/processing/Flenses';

% fill these parameters
pixsize = 0.63 ;% um %% PCO = 0.615 um  Frelon = 1.39

dist = 905;             %'nfy';% millimeters distance (sample or membrane)->whichever is downstream - detector 
Energy = 17;%12.398;            % energy in keV/mntdirect/_data_bm05_inhouse/thomas/170213lenstest1/processing/stack8
delta = 1.8777971*1e-06;%Al 17 keV = 1.8777971E-06;2.21626647.*1e-6;%Be at 12.4 keV = 2.21626647.*1e-6
maskSize = 860;   % diameter size in um you want to print the figure Vertically for the mask
rotim = -1;         % ESRF is rotated -1
savethestuff = 0;       %  boolean to save or not the stuffs
selectManuallyMask = 0;

calcRange = 1;%'all';%(1:1:12);    % range of folder you want to process 
% usefull - for the records
SetupGeometry = 'MembraneUpstream';% 'MembraneUpstream' 'SampleUpstream'
radiSource = 40200;%; in mm distance from source to fist element (membrane or sample)
lobjectdist = 400;% distance between membrane and sample in mm
%distortion loading
distorDet =  [];%open('/data/bm05/inhouse/thomas/171029lenses/processing/distortion2/disto_det.mat');
%disto = distorDet.disto;disto.ROI = distorDet.ROI;
%pixsize  = (distorDet.pixsizeV + distorDet.pixsizeH)./2;
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       calculate tau and tau2                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tau = (dist+radiSource+lobjectdist)/radiSource;
%calculate tau
if strcmp(SetupGeometry,'SampleUpstream'), tau2 = 1 + lobjectdist/radiSource;
else                                       tau2 = 1;
end;
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    Summary file                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% forget these ones

lambda = 1.2398/Energy.*1e-9;
kw2 =2*pi/lambda;

subFolders = dir([rootFolder '/*.mat']);

sampPixrate = 15;

if strcmp(calcRange,'all'), calcRange = (1:1:length(subFolders));end;
diamet = zeros(1,length(round(maskSize./pixsize./2) :-sampPixrate:sampPixrate*6));
ter   = zeros(length(calcRange),length(round(maskSize./pixsize./2) :-sampPixrate:sampPixrate*6));
verg  = zeros(length(calcRange),length(round(maskSize./pixsize./2) :-sampPixrate:sampPixrate*6));
horig = zeros(length(calcRange),length(round(maskSize./pixsize./2) :-sampPixrate:sampPixrate*6));
focV  = zeros(length(calcRange),length(round(maskSize./pixsize./2) :-sampPixrate:sampPixrate*6));  
focH  = zeros(length(calcRange),length(round(maskSize./pixsize./2) :-sampPixrate:sampPixrate*6)); 
focM  = zeros(length(calcRange),length(round(maskSize./pixsize./2) :-sampPixrate:sampPixrate*6));

if savethestuff
    if ~exist(fullfile(exportFolder,'/figures/'),'dir'),  mkdir(exportFolder,'/figures/');end;
    if ~exist(fullfile(exportFolder,'/stats/'),'dir'),    mkdir(exportFolder,'/stats/');  end;
    writeFolder = fullfile(exportFolder ,'/figures/');
    writeFolderSt = fullfile(exportFolder ,'/stats/');
end;


for kname = calcRange,
disp(   subFolders(kname).name ) ;
if savethestuff,
    exportFileSummary = fullfile(writeFolder,'parameters_2step.txt');

    dlmwrite(exportFileSummary, 'Processing parameters for second step','delimiter', '' );
    edit pmedf_readdlmwrite(exportFileSummary, ['pixsize = ' num2str(pixsize)], '-append','delimiter', '' );
    dlmwrite(exportFileSummary, ['dist = ' num2str(dist)], '-append','delimiter', '' );
    dlmwrite(exportFileSummary, ['SetupGeometry = ' num2str(SetupGeometry)], '-append','delimiter', '' );
    dlmwrite(exportFileSummary, ['radiSource = ' num2str(radiSource)], '-append','delimiter', '' );
    dlmwrite(exportFileSummary, ['lobjectdist = ' num2str(lobjectdist)], '-append','delimiter', '' );
    dlmwrite(exportFileSummary, ['Energy = ' num2str(Energy)], '-append','delimiter', '' );
    dlmwrite(exportFileSummary, ['delta = ' num2str(delta)], '-append','delimiter', '' );
    dlmwrite(exportFileSummary, ['maskSize = ' num2str(maskSize) ' um in diameter'], '-append','delimiter','');

    dlmwrite(exportFileSummary, 'Processed:', '-append','delimiter', '' );
end;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    load the 1 step processing file                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load(fullfile(rootFolder,subFolders(kname).name));

if ~isempty(distorDet),
    for k = 1:3, thephase(:,:,k) = undistorImg(thephase(:,:,k),disto); end;
end;
if rotim ~= 0, 
    thephase2 = thephase;clear thephase;
    for k = 1:4, thephase(:,:,k) = thephase2(:,:,k)';end;%flipud(rot90(thephase2(:,:,k),rotim)); end;
end;

% differential phase
thephase(:,:,1) = myerosion(thephase(:,:,1),3,3);
thephase(:,:,2) = myerosion(thephase(:,:,2),3,3);
dphaseX = thephase(:,:,1).*pixsize.*1e-3 ./dist;%vertical wft gradient % here we have radians
dphaseY = thephase(:,:,2).*pixsize.*1e-3 ./dist;%horizontal
psb = dphaseX;      psb(:,:,2) = dphaseY;


% now we have the phase gradient in radians 
[m,n,r1] = size(thephase);
mm = (1:m).*pixsize./1000;      nn = (1:n).*pixsize./1000;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           make a mask and find the center               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BUILD A MASK
[X1,Y1] = meshgrid(1:n,1:m);
X = X1 - round(n/2);  Y = Y1-round(m/2);
mask0 = (X.^2+Y.^2).^(1/2);
exponent = exp(-(((X).^2 + (Y).^2)./(2*(maskSize.*1./pixsize./2)^2)));
exponent(mask0 < (maskSize./pixsize./2)) = 1;
mask3 = (mask0 < (maskSize./pixsize./2.*1.3)).*exponent;
% Center = [159 165];% please define manually
% maskthreshold = 0.75;%0.35;% 0.5/0.45 is good unless very eak transmission inside the less
% find most intense points (proportionnaly to the mask size)
a = sort(reshape(thephase(:,:,3),1,[]),1,'descend');
maskthreshold = 1- sum(mask3(:) ==1)./numel(mask3);
a = unique(a);
b = a(round(length(a)*maskthreshold));

d1 = (thephase(:,:,3)>b) ;%+ (abs(thephase(:,:,3))>1.2)+ (abs(thephase(:,:,4))>1.2);
% for k = 1:2,% iterative morphological opening
%     d1 = (imfilter(d1,ones(5)./25,'same') >= 1);
%     d1 = (imfilter(d1,ones(5)./25,'same') > 0);
% end;

%[X,Y] = meshgrid(1:size(d1,2),1:size(d1,1));
%Center(2) = round(sum(X(d1(:)))./sum(d1(:)));
%Center(1) = round(sum(Y(d1(:)))./sum(d1(:)));

% correlation to intercept the maximum of intense points with the mask
[x_offset,y_offset] = fshift(thephase(:,:,3),mask3);    
Center(2) = x_offset + n/2;   Center(1) = y_offset + m/2;
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%            calculate wft error for defined mask size                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  mask
[X1,Y1] = meshgrid(1:n,1:m);
if ~selectManuallyMask
    %  square mask
    X = X1 - Center(2);  Y = Y1-Center(1);
    mask0 = (X.^2+Y.^2).^(1/2);
    mask2 = zeros(size(X));
    mask2 = (mask0 < (maskSize./pixsize./2));
    mask2 = mask2 == 1;
else 
    figure()
    [~,rect] = imcrop(thephase(:,:,1)./(2*max(max(thephase(:,:,1))))+.5);
    close()
    mask2 = thephase(:,:,1).*0;
    mask2(rect(2):rect(2)+rect(4),rect(1):rect(1)+rect(3)) = 1;
    mask2 = mask2 == 1;
    Center = [mean(Y1(mask2)) mean(X1(mask2))];
end;  

ps = psb(min(Y1(mask2==1)):max(Y1(mask2==1)),min(X1(mask2==1)):max(X1(mask2==1)),:);
mask1 = mask2(min(Y1(mask2==1)):max(Y1(mask2==1)),min(X1(mask2==1)):max(X1(mask2==1)));

%return;
%calculate gradient error
[p1,p2] = gradient_error(ps(:,:,1),ps(:,:,1),mask1);
[p3,p4] = gradient_error(ps(:,:,2),ps(:,:,2),mask1);
    
% calculate differential focal length
raddiff1 = 1/(mean2(diff(p1,1,2)./(pixsize*1e-6)));
raddiff2 = 1/(mean2(diff(p4,1,1)./(pixsize*1e-6)));
disp(['Differential focal length with diff fitting is in meters  ' num2str(raddiff1) '  V  ' num2str(raddiff2) ' H ']);
disp(['For a mask of diameter = ' num2str(round(maskSize)) 'um']);

% calculate a rmse moyenne
dd = ((mean2(diff(p1,1,2))./2 + mean2(diff(p4,1,1))./2 ));
p11 = p1./mean2(diff(p1,1,2)).*dd;
p44 = p4./mean2(diff(p4,1,1)).*dd;
raddiff11 = 1/(mean2(diff(p11,1,2)./(pixsize*1e-6)));
raddiff22 = 1/(mean2(diff(p44,1,1)./(pixsize*1e-6)));

%extract wavefront gradient error and put in urad (with astigmatism)

xg = (ps(:,:,1)-p11).*1e6;    xg = xg - mean(xg(:));
yg = (ps(:,:,2)-p44).*1e6;    yg = yg - mean(yg(:));
% do integration to get wavefront in meter
phi = WftSolveLSChol(yg,xg,pixsize,mask1).*1e-12;
phi = phi- mean(phi(mask1(:)));% remove piston

%extract wavefront gradient error and put in urad (without astigmatism and defocus)
xg1 = (ps(:,:,1)-p2 -p1).*1e6;    xg1 = xg1 - mean(xg1(:));
yg1 = (ps(:,:,2)-p3 -p4).*1e6;    yg1 = yg1 - mean(yg1(:));
phi1 = WftSolveLSChol(yg1,xg1,pixsize,mask1).*1e-12;
phi1 = phi1- mean(phi1(mask1(:)));% remove pistonfull map phase
phiFull = WftSolveLSChol(dphaseY - mean(dphaseY(:)),dphaseX - mean(dphaseX(:)),pixsize,ones(size(dphaseX))).*1e-12;
phiFull = phiFull- mean(phiFull(mask1(:)));% remove piston

%thickness depends on delta value
thickErr = phi./delta.*1e6;% put in microns by the way
thickErr = thickErr - mean(thickErr(mask1(:)));% remove mean
%thickness error without astimatism
thickErr1 = phi1./delta.*1e6;% put in microns by the way
thickErr1 = thickErr1 - mean(thickErr1(mask1(:)));% remove mean
% backup thickness error with astigmatism
thickErr2 = mask2.*0;
thickErr2(min(Y1(mask2)):max(Y1(mask2)),min(X1(mask2)):max(X1(mask2))) = thickErr;

% find modes

[output2] = ZernikeCalc(1:15,phi, mask1);
output2 = output2.*1e9;
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            export figures                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

exportLensFig2D

if savethestuff,
exportFileName = ['Modes_Zernike_coefficient_in_nm' '.dat']; %#ok<UNRCH>
dlmwrite(fullfile(writeFolderSt, exportFileName),strrep(subFolders(kname).name,'.mat','') ,'-append','delimiter','' );
dlmwrite(fullfile(writeFolderSt, exportFileName),output2','-append', 'delimiter', '\t' );
    
set(1,'PaperOrientation','landscape');set(3,'PaperOrientation','landscape');
set(2,'PaperOrientation','landscape');set(4,'PaperOrientation','landscape');

print(1,'-dpdf',fullfile(writeFolder,['WavefrontGradientFull_' strrep(subFolders(kname).name, '.mat','')]));
print(1,'-dpsc',fullfile(writeFolder,['WavefrontGradientFull_' strrep(subFolders(kname).name, '.mat','')]));
print(1,'-dpng',fullfile(writeFolder,['WavefrontGradientFull_' strrep(subFolders(kname).name, '.mat','')]));
%saveas(1,fullfile(writeFolder,['WavefrontGradientFull_' strrep(subFolders(kname).name, '.mat','')]));



print(2,'-dpdf',fullfile(writeFolder,['WavefrontGradientError_' strrep(subFolders(kname).name, '.mat','')]));
print(2,'-dpsc',fullfile(writeFolder,['WavefrontGradientError_' strrep(subFolders(kname).name, '.mat','')]));
print(2,'-dpng',fullfile(writeFolder,['WavefrontGradientError_' strrep(subFolders(kname).name, '.mat','')]));
%saveas(2,fullfile(writeFolder,['WavefrontGradientError_' strrep(subFolders(kname).name, '.mat','')]));


print(3,'-dpdf',fullfile(writeFolder,['Masks_' strrep(subFolders(kname).name, '.mat','')]));
print(3,'-dpsc',fullfile(writeFolder,['Masks_' strrep(subFolders(kname).name, '.mat','')]));
print(3,'-dpng',fullfile(writeFolder,['Masks_' strrep(subFolders(kname).name, '.mat','')]));
%saveas(3,fullfile(writeFolder,['Masks_' strrep(subFolders(kname).name, '.mat','')]));


print(4,'-dpdf',fullfile(writeFolder,['WavefrontThicknessError_' strrep(subFolders(kname).name, '.mat','')]));
print(4,'-dpsc',fullfile(writeFolder,['WavefrontThicknessError_' strrep(subFolders(kname).name, '.mat','')]));
print(4,'-dpng',fullfile(writeFolder,['WavefrontThicknessError_' strrep(subFolders(kname).name, '.mat','')]));
%saveas(4,fullfile(writeFolder,['WavefrontThicknessError_' strrep(subFolders(kname).name, '.mat','')]));

%print(5,'-dpng',fullfile(writeFolder,['WavefrontThicknessError_3D_' strrep(subFolders(kname).name, '.mat','')]));
end;
return;
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       vary the mask to calculate some parameters        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
r = 1;
for msz = round(maskSize./pixsize./2)+sampPixrate*6 :-sampPixrate:sampPixrate*6

    mask2 = (mask0 < msz);

    ps = psb(min(Y1(mask2==1)):max(Y1(mask2==1)),min(X1(mask2==1)):max(X1(mask2==1)),:);
    mask3 = mask2(min(Y1(mask2==1)):max(Y1(mask2==1)),min(X1(mask2==1)):max(X1(mask2==1)));
  
    
%     ps = psb(min(Y1(mask2)):max(Y1(mask2)),min(X1(mask2)):max(X1(mask2)),:);
%     mask3 = mask2(min(Y1(mask2)):max(Y1(mask2)),min(X1(mask2)):max(X1(mask2)));
    % calculate gradient error
    [p1,p2] = gradient_error(ps(:,:,1),ps(:,:,1),mask3);
    [p3,p4] = gradient_error(ps(:,:,2),ps(:,:,2),mask3);
    
    % calculate differential focal length
    raddiff1 = 1/(mean2(diff(p1,1,2)./(pixsize*1e-6)));%vertical
    raddiff2 = 1/(mean2(diff(p4,1,1)./(pixsize*1e-6)));
    disp(['Differential focal length with diff fitting is in meters  ' num2str(raddiff1) '  V  ' num2str(raddiff2) ' H ']);
    disp(['For a mask of diameter = ' num2str(round(msz.*2.*pixsize)) 'um']);
    focV(kname-calcRange(1)+1,r)  = raddiff1;    focH(kname-calcRange(1)+1,r)  = raddiff2; 

    % calculate a chiade moyenne
    dd = ((mean2(diff(p1,1,2))./2 + mean2(diff(p4,1,1))./2 ));
    p11 = p1./mean2(diff(p1,1,2)).*dd;  % vertical
    p44 = p4./mean2(diff(p4,1,1)).*dd;  % horizont
    focM(kname-calcRange(1)+1,r) = 1/(mean2(diff(p44,1,1)./(pixsize*1e-6)));

    % store some values for the records
    diamet(r) = msz.*2.*pixsize;
%     
%     xg = (ps(:,:,1)-p2 - p1).*1e6;    xg = xg - mean(xg(:));
%     yg = (ps(:,:,2)-p3 - p4).*1e6;    yg = yg - mean(yg(:));
    
    xg = (ps(:,:,1)-p11).*1e6;    xg = xg - mean(xg(:));
    yg = (ps(:,:,2)-p44).*1e6;    yg = yg - mean(yg(:));
    
%     phit = WftSolveLSChol(yg,xg,pixsize,mask3).*1e-12./delta.*1e6;
%     phit = phit- mean(phit(mask3(:)));% remove piston
%     ter(kname-calcRange(1)+1,r) = std(phit(mask3(:)));
    
    ter2(kname-calcRange(1)+1,r) = std(thickErr2(mask2(:)));
    ter = ter2;
    verg(kname-calcRange(1)+1,r)   = std(xg(mask3(:)));
    horig(kname-calcRange(1)+1,r)  = std(yg(mask3(:)));

    r = r+1;

end;
end;
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            export result matrix                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if savethestuff,
exportFileName = ['Table_FocLengthV_LensSlot' num2str(calcRange(1)-1) 'to' num2str(calcRange(end)-1) '.dat'];
dlmwrite(fullfile(writeFolderSt, exportFileName),['Table_FocLengthV_LensSlot' num2str(calcRange(1)-1) 'to' num2str(calcRange(end)-1) ] ,'delimiter','' );
dlmwrite(fullfile(writeFolderSt, exportFileName),['Mask_Diameter  '  '  Focal_distances'] ,'-append','delimiter','' );
dlmwrite(fullfile(writeFolderSt, exportFileName),[diamet' focV'],'-append', 'delimiter', '\t' );

exportFileName = ['Table_FocLengthH_LensSlot' num2str(calcRange(1)-1) 'to' num2str(calcRange(end)-1) '.dat'];
dlmwrite(fullfile(writeFolderSt, exportFileName),['Table_FocLengthH_LensSlot' num2str(calcRange(1)-1) 'to' num2str(calcRange(end)-1) ] ,'delimiter','' );
dlmwrite(fullfile(writeFolderSt, exportFileName),['Mask_Diameter  ' '  Focal_distances'] ,'-append','delimiter','' );
dlmwrite(fullfile(writeFolderSt, exportFileName),[diamet' focH'],'-append', 'delimiter', '\t' );

exportFileName = ['Table_FocLengthMean_LensSlot' num2str(calcRange(1)-1) 'to' num2str(calcRange(end)-1) '.dat'];
dlmwrite(fullfile(writeFolderSt, exportFileName),['Table_FocLengthMean_LensSlot' num2str(calcRange(1)-1) 'to' num2str(calcRange(end)-1) ] ,'delimiter','' );
dlmwrite(fullfile(writeFolderSt, exportFileName),['Mask_Diameter '  ' Focal_distances'] ,'-append','delimiter','' );
dlmwrite(fullfile(writeFolderSt, exportFileName),[diamet' focM'],'-append', 'delimiter', '\t' );

exportFileName = ['Table_ThicknessError_LensSlot' num2str(calcRange(1)-1) 'to' num2str(calcRange(end)-1) '.dat'];
dlmwrite(fullfile(writeFolderSt, exportFileName),['Table_ThicknessError_LensSlot' num2str(calcRange(1)-1) 'to' num2str(calcRange(end)-1) ] ,'delimiter','' );
dlmwrite(fullfile(writeFolderSt, exportFileName),['Mask_Diameter  '  '  ThicknessError'] ,'-append','delimiter','' );
dlmwrite(fullfile(writeFolderSt, exportFileName),[diamet' ter'],'-append', 'delimiter', '\t' );

exportFileName = ['Table_WavefrontErrorV_LensSlot' num2str(calcRange(1)-1) 'to' num2str(calcRange(end)-1) '.dat'];
dlmwrite(fullfile(writeFolderSt, exportFileName),['Table_WavefrontErrorV_LensSlot' num2str(calcRange(1)-1) 'to' num2str(calcRange(end)-1) ] ,'delimiter','' );
dlmwrite(fullfile(writeFolderSt, exportFileName),['Mask_Diameter  '  '  WavefrontErrorV'] ,'-append','delimiter','' );
dlmwrite(fullfile(writeFolderSt, exportFileName),[diamet' verg'],'-append', 'delimiter', '\t' );

exportFileName = ['Table_WavefrontErrorH_LensSlot' num2str(calcRange(1)-1) 'to' num2str(calcRange(end)-1) '.dat'];
dlmwrite(fullfile(writeFolderSt, exportFileName),['Table_WavefrontErrorH_LensSlot' num2str(calcRange(1)-1) 'to' num2str(calcRange(end)-1) ] ,'delimiter','' );
dlmwrite(fullfile(writeFolderSt, exportFileName),['Mask_Diameter  '  '  WavefrontErrorH'] ,'-append','delimiter','' );
dlmwrite(fullfile(writeFolderSt, exportFileName),[diamet' horig'],'-append', 'delimiter', '\t' );
disp(['Done with ' subFolders(kname).name]);
end;
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            print matrix graphs                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
exportMatrixFig2D

if savethestuff,
set(5,'PaperOrientation','landscape');
print('-bestfit',5,'-dpdf',fullfile(writeFolderSt,['FocalLengths-Diameter_LensSlot' num2str(calcRange(1)-1) 'to' num2str(calcRange(end)-1)]));
print('-bestfit',5,'-dpsc',fullfile(writeFolderSt,['FocalLengths-Diameter_LensSlot' num2str(calcRange(1)-1) 'to' num2str(calcRange(end)-1)]));
print(5,'-dpng',fullfile(writeFolderSt,['FocalLengths-Diameter_LensSlot' num2str(calcRange(1)-1) 'to' num2str(calcRange(end)-1)]));
saveas(5,fullfile(writeFolderSt,['FocalLengths-Diameter_LensSlot' num2str(calcRange(1)-1) 'to' num2str(calcRange(end)-1)]));
end;

