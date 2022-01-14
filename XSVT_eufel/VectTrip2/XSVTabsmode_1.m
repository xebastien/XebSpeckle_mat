% XSVT in absolute mode from two stacks of images taken at two different
% propagation distance plan

% For analysis of the EU-XFEL data
% See for the lens in the EXP hutch with no collimation upstream of the beam
% path 
%
%
% LOAD only the images
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% global varia-ble setting

close all
energy =8.7; % kev

synchrotron ='ESRF';% 'ESRF'or 'DLS'for the format of images

inputFold1 = '/data/bm05/imaging/seb/EU-XFEL/secondtrip/data2_det1/vect1';
inputFold2 = '/data/bm05/imaging/seb/EU-XFEL/secondtrip/data2_det2/vect1';

%%%%%%%%%%%%%%%%%%%5
pixsize = 1.207;            % here pixel size of det 1 since we interpolate image 2%%%1.45;%5.8;%.76;% 0.65;%6.4;%.9;%pixel size in um
maxImages = 17490;
subdiv = 12;%144*2;% chunck to use to avoid memory issues.


select_ROI = 1;% 0= full field, 1 = ROI selection, 2 = manual below 3-auto    4-reuse from 1


im_return = 0;% return afet display of the first image
filtersz = 0;% averaging filter size
myleefilter = 0;% for th1 erosion filter (morphological processing): 0 = none, value = filter threshold
undersampl = 0; % positive for binning
interpsecim = 1.25;%(1.4338673051/1.20776665922)^(1);    %          %  pixel size 2 / (pixel size 1) 1.4184
%then put a threshold number here to force the phase equal to zero in the pixel below this value
clear stack files

lambda = 12.398/energy*1e-10;
kw = 2*pi/lambda;

% -------------------------dark field and flat field images--------------%
darkfieldImagePath = [];   % [empty] for no darkfield
darkfieldImageName =  [];  %  [empty] for no file

FlatfieldImagePath = [];  % [empty] for no darkfield
FlatfieldImageName1 =  [];%'flat.edf';%'ipp10.TIF';%[empty] for no file
FlatfieldImageName2 =  [];%'flat.edf';%'ipp10.TIF';%[empty] for no file

defaultPath = '/data/bm05/imaging/seb/';% where to pick the first images

inputFoldFF = [];%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      Loading pfirst picture
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% choose a path
if (~exist('path1','var') || ~ischar(path1) || ~exist(path1,'dir')) && isempty(inputFold1) 
           path2 = uigetdir(defaultPath);
elseif isempty(inputFold1) || ~exist(inputFold1,'dir')  ,      path2 = uigetdir(path1) ;
else path2 = inputFold1;
end; %#ok<*NOSEL>

if ~ischar(path2), return;else path1 = path2;end;

display(['Path: ' path2]);fprintf('\n');
  
fprintf('\n Starting first picture \n_____');
if strcmp(synchrotron , 'ESRF')   
    rotim = -1;% rot90  the pictures
    names = [path2 '/*.edf'];
    name1 = dir(path2);
    [header1,file1] = pmedf_read(fullfile(path2,name1(3).name));
elseif strcmp(synchrotron , 'DLS')
    rotim = 0;% rot90  the pictures
    setenv('PATH2', path2)
    [~,file] = system('ls -Rp -1 $PATH2/*/**.TIF $PATH2/*.tif $PATH2/**.TIF 2>/dev/null');  
    file = cellstr(strread(file,'%s'));
    file1 = imread(file{1});
    file1 = file1(10:end,:);
end;
if myleefilter > 0, file1 = myerosion(file1, myleefilter); end; 
disp('  DONE');    fprintf('\n\n');


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
        if rotim ~= 0, darkim = rot90(darkim,rotim); end;
    end;
    if ~isempty(FlatfieldImagePath) && exist(FlatfieldImagePath,'dir')
        [flatheader,flatim1] = pmedf_read(fullfile(FlatfieldImagePath,FlatfieldImageName1));
        if rotim ~= 0, flatim1 = rot90(flatim1,rotim); end;
        if myleefilter > 0, flatim1 = myerosion(flatim1, myleefilter); end;       
        if ~isempty(darkim), flatim1 = flatim1- darkim;end; 
    end;
elseif strcmp(synchrotron , 'DLS')
    [files] = open_seq_dia(path2,ROI1);
    %-------------flat and darkfield image corection------------%
    flatim1 = []; darkim = [];
    if ~isempty(darkfieldImagePath) && exist(darkfieldImagePath,'dir')
        darkim = imread(fullfile(darkfieldImagePath,darkfieldImageName),'tif') + 1;
        if rotim ~= 0, darkim = rot90(darkim,rotim); end;
    end;
    if ~isempty(FlatfieldImagePath) && exist(FlatfieldImagePath,'dir')
        flatim1 = imread(fullfile(FlatfieldImagePath,FlatfieldImageName1),'tif');
        if rotim ~= 0, flatim1 = rot90(flatim1,rotim); end;
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
    files(k).data = files(k).data;
    if filtersz > 0 ,files(k).data = uint16(imfilter(files(k).data,ones(filtersz)./(filtersz).^2,'replicate','same'));end;
    if exist('rotim','var')   && (rotim ~= 0), files(k).data = rot90(files(k).data,rotim );end;
    if ~isempty(darkim),  files(k).data = single(files(k).data) - darkim;end;
    if ~isempty(flatim1), files(k).data = single(files(k).data) ./ single(flatim1);end;
    if undersampl > 0,files(k).data = files(k).data(1:undersampl:end,1:undersampl:end); end;
    progmeter(k/nImages);
end;
if undersampl > 0, pixsize = pixsize *undersampl;end;

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
%                  det2 image loading
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
if isempty(inputFold2) || ~exist(inputFold2,'dir') ,      
    path2 = uigetdir(fullfile(path1,'/..'));
else 
    path2 = inputFold2;
end;
disp(path2);
filesSample = files;% save for the 
clear files
% open the files from the path
fprintf('\n Starting loading Flatfield pictures\n_____');
if strcmp(synchrotron , 'ESRF'),    
    [files] = open_seq(path2);
    %-------------flat and darkfield image corection------------%
    flatim2 = []; darkim = [];
    if ~isempty(darkfieldImagePath) && exist(darkfieldImagePath,'dir')
        [darkheader, darkim] = pmedf_read(fullfile(darkfieldImagePath,darkfieldImageName));
        if rotim ~= 0, darkim = rot90(darkim,rotim); end;
    end;
    if ~isempty(FlatfieldImagePath) && exist(FlatfieldImagePath,'dir')
        [flatheader,flatim2] = pmedf_read(fullfile(FlatfieldImagePath,FlatfieldImageName2));
        if rotim ~= 0, flatim1 = rot90(flatim2,rotim); end;
        if myleefilter > 0, flatim2 = myerosion(flatim2, myleefilter); end;       
        if ~isempty(darkim), flatim2 = flatim2- darkim;end; 
    end;
elseif strcmp(synchrotron , 'DLS'),
   [files] = open_seq_dia(path2);
   %-------------flat and darkfield image corection------------%
   if ~isempty(FlatfieldImagePath) && exist(FlatfieldImagePath,'dir')
        flatim2 = imread(fullfile(FlatfieldImagePath,FlatfieldImageName2),'tif');
        if rotim > 0, flatim2 = rot90(flatim2,rotim); end;
        if myleefilter > 0, flatim2 = myerosion(flatim2, myleefilter); end;       
        if ~isempty(darkim), flatim2 = flatim2- darkim;end; 
   end;
end;
disp(' Flatfield images loaded');

% error checking
if ~isstruct(files) || isempty(files(1)), error('No FlatField Spictures in memory: path probably not correct');end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[m1 ,n1] = size(files(1).data); 
[XI,YI] = meshgrid(linspace(1,n1,round(n1/interpsecim)),linspace(1,m1,round(m1/interpsecim)));




%%% reduce to ROI and process
for k = 1 : 1 : nImages, 
    files(k).data = files(k).data;
    if myleefilter > 0, files(k).data = myerosion(files(k).data, myleefilter,'silent'); end; 
    if filtersz > 0 ,files(k).data = uint16(imfilter(files(k).data,ones(filtersz)./(filtersz).^2,'replicate','same'));end;
    if exist('rotim','var')   && (rotim ~= 0), files(k).data = rot90(files(k).data,rotim );end;
    if ~isempty(darkim), files(k).data = single(files(k).data) - darkim;end;
    if ~isempty(flatim1), files(k).data = single(files(k).data) ./ single(flatim1);end;
    if undersampl > 0,files(k).data = files(k).data(1:undersampl:end,1:undersampl:end); end;
    if interpsecim > 0, files(k).data = interp2(single(files(k).data),XI,YI,'linear'); end
    progmeter(k/nImages);
end;

figure(6)
subplot(1,2,2)
imagesc(files(1).data)
colorbar
colormap(gray)
title('First image ref')

%% %%%%%%%%%%%%%%%%%%%%%%%reorgnize name and size%%%%%%%%%%%%%%%%%%%%%%%%%
files1 = filesSample;
files2 = files;

[m1,n1] = size(files1(1).data);
[m2,n2] = size(files2(1).data);
clear files filesSample

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% reshape the data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
zonebounds1 = round(linspace(1,m1+1,subdiv+1));
zonebounds2 = round(linspace(1,m2+1,subdiv+1));
piece_stack_sample = cell(1,subdiv);
piece_stack_ref = cell(1,subdiv);
for k = 1 :subdiv
    zone1 = zonebounds1(k):1:zonebounds1(k+1)-1;
    zone2 = zonebounds2(k):1:zonebounds2(k+1)-1;
    piece_stack_sample{k} = zeros(length(zone1) , n1, nImages );
    piece_stack_ref{k} = zeros(length(zone2) , n2, nImages );
    for pp = 1 : nImages, 
        imgs = files1(pp).data;
        piece_stack_sample{k}(:,:,pp) = imgs(zone1,:);
        imgr = files2(pp).data;
        piece_stack_ref{k}(:,:,pp) = imgr(zone2,:);
        piece_stack_sample{k} = uint16(piece_stack_sample{k});
        piece_stack_ref{k} = uint16(piece_stack_ref{k});
    end;
end;


stsample = cat(1,piece_stack_sample{:}); 
stref = cat(1,piece_stack_ref{:});
%clear piece_stack_sample piece_stack_ref
[m,n,r] = size(stsample);
file1 = mean(stsample,3);
return;