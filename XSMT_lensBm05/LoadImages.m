%% Load two structures of images filesSample and files (ref ones)
%   by S. Berujon April 2016
%   contact: dont@contact.me
%   this script only load the images in memory
%   
%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% global varia-ble setting
synchrotron ='ESRF';% 'ESRF'or 'DLS'

%%%%%%%%%%%%%%1%%%%%ROi selection%%%%%%%%%%%%%%%%%
select_ROI = 0;% 0= full field, 1 = ROI selection, 2 = manual below 3-auto    4-reuse from 1
ROI1 = [1 900 1 950];                             % 
%%%%%%%%%%filters%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rotim = 0;      % rot90  the pictures
filtersz = 0;   % averaging filter size
myleefilter = 0;% for the erosion filter (morphological processing): 0 = none, value = filter threshold
undersampl = 0; % positive for binning
%then put a threshold number here to force the phase equal to zero in the pixel below this value
% PCO = 475 depending on exposure, MFDI = 105
h = fspecial('gaussian', 3);% 0 for no filtering

im_return = 0;% return afet display of the first image
maxImages = 5555;% number of max images - ignore images beyond this number
subdiv = 12; % number chuncks to create to avoid memory issues and compute in parallel 
%% =============================== flat and dark ===============
% -------------------------dark field and flat field images--------------%
darkfieldImagePath = [];%'/data/visitor/mi1227/id16b/processing/DarkFlat/';%';% [empty] for no darkfield
darkfieldImageName =  [];%'darkF.edf';%'ipp1.TIF';%[empty] for no file

FlatfieldImagePath = [];%'/data/id06/inhouse/2016/run2/160405_nfs/processing/Lens/';%'/data/id06/inhouse/2016/run2/160405_nfs/processing/Lens/';%'/data/visitor/mi1227/id16b/processing';% [empty] for no darkfield
FlatfieldImageName1 = [];% 'flatAc18.edf';%'flat_dist=300.edf';%'flat.edf';%'ipp10.TIF';%[empty] for no file
FlatfieldImageName2 =  [];%'flatAc18.edf';%'flat_dist=300.edf';%'flat.edf';%'ipp10.TIF';%[empty] for no file

defaultPath = '/data/visitor/in994';% where to pick the first images

%% =============================== distortion ===============
distoFileName = [];%'/users/berujon/figures/XST_SG_ID06/disto/disto_det1.mat';

% only used if called from another script to avoid human clicking interaction
inputFold1 = [];
inputFoldFF = [];
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      Loading pfirst picture
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%close all;
clear stack files;
% choose a path
if (~exist('path1','var') || ~ischar(path1) || ~exist(path1,'dir')) && isempty(inputFold1) 
           path2 = uigetdir(defaultPath);
elseif isempty(inputFold1) || ~exist(inputFold1,'dir')  ,      path2 = uigetdir(path1) ;
else path2 = inputFold1;
end;

if ~ischar(path2), return;else path1 = path2;end;


display(['Path: ' path2]);fprintf('\n');
  
fprintf('\n Starting first picture \n_____');
if strcmp(synchrotron , 'ESRF')   
    names = [path2 '/*.edf'];
    name1 = dir(path2);
    [header1,file1] = pmedf_read(fullfile(path2,name1(6).name));
elseif strcmp(synchrotron , 'DLS')
    setenv('PATH2', path2)
    [~,file] = system('ls -Rp -1 $PATH2/*/**.TIF $PATH2/*.tif $PATH2/**.TIF 2>/dev/null');  
    file = cellstr(strread(file,'%s')); %#ok<FPARK>
    file1 = imread(file{6});
    file1 = file1(10:end,:);
end;
if myleefilter > 0, file1 = myerosion(file1, myleefilter); end; 
disp('  DONE');    fprintf('\n\n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                  Load distortion files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isempty(distoFileName),
   fileStruct = load(distoFileName);
   ROIdisto = fileStruct.ROI;     disto = fileStruct.disto;
   
   [file1] = undistorImg(file1,disto,ROIdisto);
end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    selection of the ROI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% we consider here that the spot of interest doesn't move or is bigger than
% the field of view
figure(7)
[m, n] = size(file1);
if select_ROI == 0
    ROI1 = ones(1,4);
    ROI1(2) = size(file1,1);
    ROI1(4) = size(file1,2);
elseif select_ROI == 1
    ROI = zeros(1,4);% dim1down -> dim1 up , dim2 down -> dim2 up
    [A,rect] = imcrop(single(file1)./max(single(file1(:))));    
    ROI1 = round([rect(2) rect(2)+rect(4) rect(1) rect(1)+rect(3)]);
    ROI2 = ROI1;
elseif select_ROI == 2    
    ROI = ROI1;
elseif select_ROI == 3
    [ROI5, ROI6, ROI7] = roi_calc_v2(file1,40,0.1,[15 15 15 15], 300);
    ROI1 = ROI7;
elseif select_ROI == 4 
    if ~exist('ROI2','var');display('No previous Roi in memory - returning');return;end;
    ROI1 = ROI2;%unwrapped2  = unwrap(unwrap(phaseMap,[],2 ));
end;
disp(['ROI use is   ' num2str(ROI1)]);
disp('------')

close(7)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    open the files of the sample
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\n Starting loading pictures\n_____');
if strcmp(synchrotron , 'ESRF')    
    clear files
    [files] = open_seq(path1);
    %-------------flat and darkfield image corection------------%
    flatim1 = []; darkim = [];
    if ~isempty(darkfieldImagePath) && exist(darkfieldImagePath,'dir')
        [darkheader, darkim] = pmedf_read(fullfile(darkfieldImagePath,darkfieldImageName));
        % darkim = darkim(ROI1(1):ROI1(2),ROI1(3):ROI1(4));
        % if rotim ~= 0, darkim = rot90(darkim,rotim); end;
    end;
    if ~isempty(FlatfieldImagePath) && exist(FlatfieldImagePath,'dir')
        [flatheader,flatim1] = pmedf_read(fullfile(FlatfieldImagePath,FlatfieldImageName1));
        % flatim1 = flatim1(ROI1(1):ROI1(2),ROI1(3):ROI1(4));
        % if rotim ~= 0, flatim1 = rot90(flatim1,rotim); end;
        if myleefilter > 0, flatim1 = myerosion(flatim1, myleefilter); end;       
        if ~isempty(darkim), flatim1 = flatim1- darkim;end; 
    end;
elseif strcmp(synchrotron , 'DLS')
    [files] = open_seq_dia(path2);
    %-------------flat and darkfield image corection------------%
    flatim1 = []; darkim = [];
    if ~isempty(darkfieldImagePath) && exist(darkfieldImagePath,'dir')
        darkim = imread(fullfile(darkfieldImagePath,darkfieldImageName),'tif') + 1;
        % darkim = darkim(ROI1(1):ROI1(2),ROI1(3):ROI1(4));
        % if rotim ~= 0, darkim = rot90(darkim,rotim); end;
    end;
    if ~isempty(FlatfieldImagePath) && exist(FlatfieldImagePath,'dir')
        flatim1 = imread(fullfile(FlatfieldImagePath,FlatfieldImageName1),'tif');
        % flatim1 = flatim1(ROI1(1):ROI1(2),ROI1(3):ROI1(4));
        % if rotim ~= 0, flatim1 = rot90(flatim1,rotim); end;
        if myleefilter > 0, flatim1 = myerosion(flatim1, myleefilter); end;       
        if ~isempty(darkim), flatim1 = flatim1- darkim;end; 
    end;
end;
disp('  DONE');

% error checking
if ~isstruct(files) || isempty(files(1)), error('No 1pictures in memory: path probably not correct');end;
nImages = length(files);
if nImages > maxImages, nImages = maxImages;end;
disp(['Number of pictures of the sample = ' num2str(nImages) ' or ' num2str(sqrt(nImages)) ' ^2'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% reduce to ROI and process
for k = 1 : 1 : nImages, 
    %files(k).data = files(k).data;
    if myleefilter > 0, files(k).data = myerosion(files(k).data, myleefilter,'silent'); end; 
    if ~isempty(darkim), files(k).data = single(files(k).data) - darkim;end;
    if ~isempty(flatim1), files(k).data = single(files(k).data) ./ single(flatim1);end;
    if ~isempty(distoFileName),files(k).data = undistorImg(files(k).data,disto,ROIdisto); end;
    files(k).data = files(k).data(ROI1(1):ROI1(2),ROI1(3):ROI1(4));
    if exist('rotim','var')   && (rotim ~= 0), files(k).data = rot90(files(k).data,rotim );end;
    if filtersz > 0 ,files(k).data = uint16(imfilter(files(k).data,ones(filtersz)./(filtersz).^2,'replicate','same'));end;
    if undersampl > 0,imfiltered = imfilter(files(k).data,ones(3)./9,'same'); files(k).data= imfiltered(1:undersampl:end,1:undersampl:end); end;
    progmeter(k/nImages);
%     d = fft2(files(k).data);
%     d(m1-f1f:m1+f1f,n1-f2f:n1+f2f) = 0;
%     files(k).data = ifft2(d,'symmetric');
    if h(1) ~= 0, files(k).data = imfilter(files(k).data,h,'same');    end;
end;


[m1, n1] = size(files(1).data);
% _______________________Plot the first image______________________________%
disp('***************');
[m, n] = size(files(1).data);
figure(6)
subplot(1,2,1)
imagesc(rot90(file1,rotim))
colorbar
colormap(gray)
title('First image')
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                  reference image loading
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
if isempty(inputFold1) || ~exist(inputFold1,'dir') ,      
    path2 = uigetdir(fullfile(path1,'/..'));
else 
    path2 = inputFoldFF;
end;
disp(path2);
filesSample = files;% save for the 
clear files
% open the files from the path
fprintf('\n Starting loading Flatfield pictures\n_____');
if strcmp(synchrotron , 'ESRF'),    [files] = open_seq(path2);
elseif strcmp(synchrotron , 'DLS'),
   [files] = open_seq_dia(path2);
   %-------------flat and darkfield image corection------------%
   if ~isempty(FlatfieldImagePath) && exist(FlatfieldImagePath,'dir')
        flatim2 = imread(fullfile(FlatfieldImagePath,FlatfieldImageName2),'tif');
        % flatim2 = flatim2(ROI1(1):ROI1(2),ROI1(3):ROI1(4));
        % if rotim > 0, flatim2 = rot90(flatim2,rotim); end;
        if myleefilter > 0, flatim2 = myerosion(flatim2, myleefilter); end;       
        if ~isempty(darkim), flatim2 = flatim2- darkim;end; 
   end;
end;
disp(' Flatfield images loaded');

% error checking
if ~isstruct(files) || isempty(files(1)), error('No FlatField Spictures in memory: path probably not correct');end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% reduce to ROI and process
for k = 1 : 1 : length(files), 
    if myleefilter > 0, files(k).data = myerosion(files(k).data, myleefilter,'silent'); end; 
    if ~isempty(darkim), files(k).data = single(files(k).data) - darkim;end;
    if ~isempty(flatim1), files(k).data = single(files(k).data) ./ single(flatim1);end;
    if ~isempty(distoFileName),files(k).data = undistorImg(files(k).data,disto,ROIdisto); end;
    files(k).data = files(k).data(ROI1(1):ROI1(2),ROI1(3):ROI1(4));
    if exist('rotim','var')   && (rotim ~= 0), files(k).data = rot90(files(k).data,rotim );end;
    if filtersz > 0 ,files(k).data = uint16(imfilter(files(k).data,ones(filtersz)./(filtersz).^2,'replicate','same'));end;
    %if undersampl > 0,files(k).data = files(k).data(1:undersampl:end,1:undersampl:end); end;
    if undersampl > 0,imfiltered = imfilter(files(k).data,ones(3)./9,'same'); files(k).data= imfiltered(1:undersampl:end,1:undersampl:end); end;
    progmeter(k/nImages);
%     d = fft2(files(k).data);
%     d(m1-f1f:m1+f1f,n1-f2f:n1+f2f) = 0;
%     files(k).data = ifft2(d,'symmetric');
    if h(1) ~= 0,files(k).data = imfilter(files(k).data,h,'same'); end;
end;

figure(6)
subplot(1,2,2)
imagesc(files(1).data)
colorbar
colormap(gray)
title('First image ref')
% %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% flat correct %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% avgImg = zeros(size(files(1).data));
% for k = 1:nImages, avgImg = avgImg + single(files(k).data)./nImages;    end;
% for k = 1:nImages, files(k).data = single(files(k).data)./avgImg;end;
% for k = 1:nImages, filesSample(k).data = single(filesSample(k).data)./avgImg;end;

nImagesR = length(files);
lim1 = (nImagesR - nImages)/2;


for  k = 1:length(files), 
    files(k).data = (single(files(k).data) - single(mean(files(k).data(:))))./std(single(files(k).data(:)));
    filesSample2(k).data = files(k).data.*0;
end;


for  k = 1:length(filesSample), 
    filesSample(k).data = (single(filesSample(k).data) - single(mean(filesSample(k).data(:))))./std(single(filesSample(k).data(:)));
end;


for k = (lim1+1:1:lim1+nImages)
    filesSample2(k).data = filesSample(k-lim1).data;
end;

filesSample = filesSample2;
clear filesSample2;
nImages = nImagesR;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% reshape the data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
zonebounds = round(linspace(1,m1+1,subdiv+1));
piece_stack_sample = cell(1,subdiv);
piece_stack_ref = cell(1,subdiv);
for k = 1 :subdiv
    zone = zonebounds(k):1:zonebounds(k+1)-1;
    piece_stack_sample{k} = zeros(length(zone) , n1, nImages );
    piece_stack_ref{k} = zeros(length(zone) , n1, nImages );
    for pp = 1 : nImages, 
        imgs = filesSample(pp).data;
        piece_stack_sample{k}(:,:,pp) = imgs(zone,:);
        imgr = filesSample(pp).data.*0;
        piece_stack_ref{k}(:,:,pp) = imgr(zone,:);
    end;
end;


return;