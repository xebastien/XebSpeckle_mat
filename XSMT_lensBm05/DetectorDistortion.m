
% Calculation of the pixel size from a scan of the detector position % 
% + CALULCATION of the detector distorio
% code started 12/03/2014
% INPUT : at top of the script
%           pick a folder with the scan images -> mesh scan (20 um)step of
%           the detector position
%           
% fast axis is vertical in the referencial of the image(horizontal on
% frelon)
% Roi of interest = full
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
synchrotron ='ESRF';% 'ESRF'or 'DLS' 
%---------------
grid_resol = 2;% resolution of the grid where the speckle moves will be calculated
%%%%%%%%%%% physical parameters %%%6%%%%%%%%%%%%
corrsize = 13;% shalf the size of the cross-correlation subsets
motormneV = 'y3';% 'diff2z';% mnemonique of the motor scanned (before image rotation in matlab )
motormneH = 'frelz';%'diff2y';% mnemonique of the motor scanned
%orientationscan = 'horizontal';% 'vertical'or 'horizontal' => you will have to run this script twice for complete characterization distortion
kindofAvg = 'fancy';% or'median' or average or fancy
id22 = 0;% if data are coming from id22, a line is added to compensate for the detector line missing
intercorrStep = 1;% how far on the mesh do you want to correlate the various pictures
%---------------for ESRF x3 use -----------------------------------%
ROImode = 'input';%'manual'(drag and drop),'repeat' or 'input' 
%[top down left right]  = [mindim1 maxdim1 mindifi
ROIinput = [31 (2048-30) 31 (2048-30)];
% ROIinput = [61        1988          61        1988];
% ROIinput = [31        1988          61        1988];
%------------------------instructions to process--------------------------%
stop_imshow = 0;% one if you want to stop the function after having dispayer the first images
myleefilter = 0 ; % activate my noise remover filter (morphological erosion)
noreload = 0;% do ypu want to use the images already in memory
filtersz = 0 ;% averaging filter
undersampl = 0;% binning factor
% -------------------------dark field and flat field images--------------%
darkfieldImagePath = [];% [empty] for no darkfield
darkfieldImageName = [];% [empty] for no file

% FlatfieldImagePath = [];%'/data/bm05/inhouse/thomas/171211lenses/processing/detdisto';%xst/xst_meshdet1/';% [empty] for no darkfield
% FlatfieldImageName1 =  'avim';%'distoFlat.edf';%'ipp10.TIF';%[empty] for no file
% FlatfieldImageName1 =  'avim' to use the median image 
%savePlaceName = '/data/bm05/inhouse/thomas/171211lenses/processing/detdisto/GradDist_avg.mat';

defaultPath = '/data/bm05/inhouse/';
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      Loading the first picture
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('files','var'), noreload = 0;end;
if ~noreload
% choose a path
% if ~exist('path1','var') || ~ischar(path1) || ~exist(path1,'dir')
%            path2 = uigetdir(defaultPath);
% else       path2 = uigetdir(path1);
% end;
intercorrStep2 = intercorrStep;
if ~ischar(path2), return;else path1 = path2;end;
display(['Path: ' path2]);
  fprintf('\n');
% Load the first image
fprintf('\n Starting first picture \n_____');
if strcmp(synchrotron , 'ESRF')    
    names = [path2 '/*.edf'];
    name1 = dir(names);
    [header1, file1] = pmedf_read(fullfile(path2,name1(3).name));
        %-------------flat and darkfield image corection------------%
    flatim1 = ones(size(file1)); darkim = zeros(size(file1));
    if ~isempty(darkfieldImagePath) && exist(darkfieldImagePath,'dir')
        [darkheader, darkim] = pmedf_read(fullfile(darkfieldImagePath,darkfieldImageName));
    end;
    if ~isempty(FlatfieldImagePath) && exist(FlatfieldImagePath,'dir') &&  ~strcmp(FlatfieldImageName1,'avim')
        [flatheader,flatim1] = pmedf_read(fullfile(FlatfieldImagePath,FlatfieldImageName1));
        if myleefilter > 0, flatim1 = myerosion(flatim1, myleefilter); end;       
        if ~isempty(darkim),flatim1 = flatim1- darkim; end; 
    end;
elseif strcmp(synchrotron , 'DLS')
    names = [path2 '/ipp*'];
    name1 = dir(path2);
    file1 = imread(fullfile(path2,name1(3).name));
end;
disp('  DONE');    fprintf('\n\n');
file1 = (single(file1)-single(darkim))./(single(flatim1));

if undersampl > 0,file1 = file1(1:undersampl:end,1:undersampl:end); end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    selection of the ROI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% we consider here that the spot of interest doesn't move or is bigger than
% the field of view
if ~exist('ROIbkp','var'), ROIbkp = []; end; 
[ROI] = ROIselect(file1,ROImode,ROIinput,ROIbkp);
ROIbkp = ROI;% save for rerun
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    load all the images from the path (must have a square
%                    number)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\n Starting loading pictures\n_____');
if strcmp(synchrotron , 'ESRF')    
    [files] = open_seq(path1);
elseif strcmp(synchrotron , 'DLS')
    [files] = open_seq_dia(path2);
end;
disp('  DONE');

% error checking
if ~isstruct(files) || isempty(files(1)), error('No 1pictures in memory: path probably not correct');end;

nImages = length(files);
disp(['Number of pictures = ' num2str(nImages)])
disp('     ');
% display the first image of the stack
[m, n] = size(files(1).data);


% calculate average flat image if option acticated
if strcmp(FlatfieldImageName1,'avim'), 
    flatim1 = single(files(1).data).*0;
    for k = 1:nImages, flatim1 = flatim1 + single(files(k).data)./nImages;end;
end;
%for k = 1:nImages, files(k).data = undistorImg(single(files(k).data),disto);end;

figure(6)
subplot(1,2,1)
imagesc(file1)
colorbar
colormap(gray)
subplot(1,2,2)
imagesc(flatim1)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                     reduce to ROI and filter
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xdet = zeros(1,nImages);% store x detector position of the images
ydet = zeros(1,nImages);% store y detector position of the images
disto.ROI = ROI;
for k = 1 :nImages, 
    if myleefilter, files(k).data = myerosion(files(k).data, myleefilter,'silent');  end;
    if ~isempty(darkim),  files(k).data = single(files(k).data) - single(darkim);  end;
    if ~isempty(flatim1), files(k).data = single(files(k).data) ./ single(flatim1);end;%flat already dark corrected
    if id22, % correct for chip little soucis
        files(k).data = [files(k).data(:,1:1024) files(k).data(:,1024) files(k).data(:,1025:end)];
    end;
    if filtersz > 0 ,files(k).data = imfilter(files(k).data,ones(filtersz)./(filtersz).^2,'replicate','same');end;
    if undersampl > 0,files(k).data = files(k).data(1:undersampl:end,1:undersampl:end); end;
    
    % [files(k).data] = undistorImg(files(k).data,disto);
    progmeter(k/nImages);
    xdet(k) = pmedf_findPos( files(k).header, 'motor', motormneH );
    ydet(k) = pmedf_findPos( files(k).header, 'motor', motormneV );
end;


if mode(diff(xdet)) == 0,
    fastAxis = 'YY';
else 
    fastAxis = 'XX';
end;



end;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     apply xst to recover dtector distortion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
display('Calculating speckle moves');
% define the images for which we are going to calculate the move (exclude
% the images at the edge of the mesh)
rtnImages = sqrt(nImages);
if (strcmp(orientationscan,'horizontal') && strcmp(fastAxis,'XX')) || (strcmp(orientationscan,'vertical') && strcmp(fastAxis,'YY'))
ksel = (1:nImages);
    ksel = ksel(mod(ksel,rtnImages)~=0 );% remove the images on the 
    if intercorrStep >1
    for kk = 1:1:(intercorrStep-1)
        ksel = ksel(mod(ksel,rtnImages)~=(rtnImages-kk) );
    end;
    end;
else
    intercorrStep = intercorrStep*rtnImages;
    ksel = 1 : 1 : nImages -intercorrStep;% remove n last lines
end;

SpeckDis1 = cell(1,length(ksel));
SpeckDis2 = cell(1,length(ksel));
% we now calculate in the teo directions for later averaging
figure(8)
for k = (1:1:length(ksel))%
    tic
    [SpeckDis2{k},~] = cross_spot4(files(ksel(k)+intercorrStep).data,files(ksel(k)).data,ROI,ROI,grid_resol,2, ones(size(file1)),corrsize); 
    toc
    subplot(1,2,1)
    imagesc(SpeckDis2{k}(:,:,1))
    title(['Avg = ' num2str(mean(reshape(SpeckDis2{k}(:,:,1),1,[]))) '  std = ' num2str(std(reshape(SpeckDis2{k}(:,:,1),1,[])))])
    subplot(1,2,2)
    imagesc(SpeckDis2{k}(:,:,2))
    title(['Avg = ' num2str(mean(reshape(SpeckDis2{k}(:,:,2),1,[]))) '  std = ' num2str(std(reshape(SpeckDis2{k}(:,:,2),1,[])))])
    drawnow
    disp(['Calc '  num2str(k) ' of ' num2str(length(ksel)) ' part 1']);
end;

for k = 1:1:length(ksel)%
    tic
    [SpeckDis1{k},~] = cross_spot4(files(ksel(k)).data,files(ksel(k)+intercorrStep).data,ROI,ROI,grid_resol,2, ones(size(file1)),corrsize); 
    toc
    subplot(1,2,1)
    imagesc(SpeckDis1{k}(:,:,1))
    title(['Avg = ' num2str(mean(reshape(SpeckDis1{k}(:,:,1),1,[]))) '  std = ' num2str(std(reshape(SpeckDis1{k}(:,:,1),1,[])))])
    subplot(1,2,2)
    imagesc(SpeckDis1{k}(:,:,2))
    title(['Avg = ' num2str(mean(reshape(SpeckDis1{k}(:,:,2),1,[]))) '  std = ' num2str(std(reshape(SpeckDis1{k}(:,:,2),1,[])))])
    drawnow
    disp(['Calc '  num2str(k) ' of ' num2str(length(ksel)) ' part 2']);
end;
% prepare variable for collecting cell
SpeckDisV = zeros(size(SpeckDis1{1}(:,:,1)));
SpeckDisH = zeros(size(SpeckDis1{1}(:,:,1)));
SpeckDisV1 = zeros(size(SpeckDis1{1}(:,:,1)));
SpeckDisH1 = zeros(size(SpeckDis1{1}(:,:,1)));
% build sstacks
for k = 1:1:length(SpeckDis1)
    SpeckDisV1 = SpeckDisV1 + (SpeckDis1{k}(:,:,1)-SpeckDis2{k}(:,:,1))./2;
    SpeckDisH1 = SpeckDisH1 + (SpeckDis1{k}(:,:,2)-SpeckDis2{k}(:,:,2))./2;
    
    SpeckDisV(:,:,k) = (SpeckDis1{k}(:,:,1)-SpeckDis2{k}(:,:,1))./2;
    SpeckDisH(:,:,k) = (SpeckDis1{k}(:,:,2)-SpeckDis2{k}(:,:,2))./2;
end;

disp(['Number of correlation images  : ' num2str(length(SpeckDis1))]);
% calculate median value
SpeckDisV = sort(SpeckDisV,3);
SpeckDisH = sort(SpeckDisH,3);
ll = round(size(SpeckDisV,3)/3);


SpeckDisV2 = mean(SpeckDisV(:,:,ll:end-ll),3);
SpeckDisH2 = mean(SpeckDisH(:,:,ll:end-ll),3);

SpeckDisV = median(SpeckDisV,3);
SpeckDisH = median(SpeckDisH,3);

%calculate average value
SpeckDisV1 = SpeckDisV1./length(SpeckDis1);%(nImages -intercorrStep);
SpeckDisH1 = SpeckDisH1./length(SpeckDis1);%(nImages -intercorrStep);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% median values
if strcmp (kindofAvg,'median')
    ssV = SpeckDisV;                ssH = SpeckDisH;
elseif strcmp (kindofAvg,'average')% or average values
    ssV = SpeckDisV1;               ssH = SpeckDisH1;
elseif strcmp (kindofAvg,'fancy')% or average values
    ssV = SpeckDisV2;               ssH = SpeckDisH2;
end;



maskI = flatim1 >= (mean(flatim1(:)/2))-1;
maskI = maskI(ROIinput(1):ROIinput(2),ROIinput(3):ROIinput(4));
maskI = maskI(1:size(ssV,1),1:size(ssV,2));

if strcmp(orientationscan,'horizontal')
    if mean(ssH(:)) < 0, ssV = -ssV;ssH = -ssH;end;
    detstep = abs(diff(xdet));detstep = median(detstep(detstep > 0));
    pixsizeH = abs(detstep.*intercorrStep2 .*1000./mean(ssH(:)));
    disp(['Horizontal pixel size = ' num2str(pixsizeH) ' um (before image rotation)']);
    sdhv = (ssV - mean(ssV(:)))./mean(ssH(:));
    sdhh = (ssH - mean(ssH(:)))./mean(ssH(:));
    if ~exist(savePlaceName,'file')
        save(savePlaceName,'sdhh','sdhv','pixsizeH','ROI');
    else
        save(savePlaceName,'sdhh','sdhv','pixsizeH','ROI','-append')
    end;
elseif strcmp(orientationscan,'vertical')
    if mean(ssV(:)) < 0, ssV = -ssV;ssH = -ssH;end;
    detstep = abs(diff(ydet));detstep = median(detstep(detstep > 0));
    pixsizeV = detstep.*intercorrStep2.*1000./abs(mean(ssV(maskI)));
    disp(['Vertical pixel size = ' num2str(pixsizeV)   ' um (before image rotation)']);
    sdvv = (ssV - mean(ssV(:)))./mean(ssV(:));
    sdvh = (ssH - mean(ssH(:)))./mean(ssV(:));
    if ~exist(savePlaceName,'file')
        save(savePlaceName,'sdvh','sdvv','pixsizeV','ROI');
    else
        save(savePlaceName,'sdvh','sdvv','pixsizeV','ROI','-append')
    end;
end;
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


save(savePlaceName,'maskI','-append')


figure(5)
subplot(1,2,1)
imagesc(ssV)
title('ssV')
colorbar

subplot(1,2,2)
imagesc(ssH)
title('ssH')
colorbar