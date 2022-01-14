%% Speckle darkfield and scatterinf imaging from membrane scan
%   by S. Berujon April 2012
%   contact: dont@contact.me
%   this script only load the images in memory
%   
%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% global varia-ble setting

close all
energy = 17;%kev

synchrotron ='ESRF';% 'ESRF'or 'DLS'

%%%%%%%%%%%%%%%%%%%5
pixsize = 1.6;%5.8;%.76;% 0.65;%6.4;%.9;%pixel size in um
maxImages = 17490;
subdiv = 12;%144*2;% chunck to use to avoid memory issues.


select_ROI = 1;% 0= full field, 1 = ROI selection, 2 = manual below 3-auto    4-reuse from 1
% ROI1 = [978 1756 455 3500];% G1G2=217
ROI1 = [1 900 1 950];% 



im_return = 0;% return afet display of the first image
filtersz = 0;%15;% averaging filter size
myleefilter = 0;% for the erosion filter (morphological processing): 0 = none, value = filter threshold
undersampl = 2; % positive for binning
%then put a threshold number here to force the phase equal to zero in the pixel below this value
% PCO = 475 depending on exposure, MFDI = 105
clear stack files
%/mntdirect/_data_visitor/mi1041/bm05/phstepping_7thorde
lambda = 12.398/energy*1e-10;
kw = 2*pi/lambda;


% -------------------------dark field and flat field images--------------%
darkfieldImagePath = [];%'/data/visitor/mi1227/id16b/processing';% [empty] for no darkfield
darkfieldImageName =  [];%'dark.edf';%'ipp1.TIF';%[empty] for no file

FlatfieldImagePath = '/data/visitor/md950/id17/cerveau_manu_speckles_radio/';%'/data/visitor/mi1227/id16b/processing';% [empty] for no darkfield
FlatfieldImageName1 =  'flat.edf';%'ipp10.TIF';%[empty] for no file
FlatfieldImageName2 =  [];%'flat.edf';%'ipp10.TIF';%[empty] for no file

defaultPath = '/dls/b16/data/2013/nt5871-2';% where to pick the first images


inputFold1 = [];%'/mntdirect/_data_bm05_inhouse/seb/140709_CohColloids/Ant_Piezo_P1200/MS';% empty to chose manually
inputFoldFF = '/mntdirect/_data_bm05_inhouse/seb/140709_CohColloids/Ant_Piezo_P1200/MF';
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      Loading pfirst picture
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    rotim = 1;% rot90  the pictures
    names = [path2 '/*.edf'];
    name1 = dir(path2);
    [header1,file1] = pmedf_read(fullfile(path2,name1(3).name));
elseif strcmp(synchrotron , 'DLS')
    rotim = 0;% rot90  the pictures
    setenv('PATH2', path2)
    [~,file] = system('ls -Rp -1 $PATH2/*/**.TIF $PATH2/*.tif $PATH2/**.TIF $PATH2/*/**.TIFF $PATH2/*.tiff $PATH2/**.TIFF 2>/dev/null');  
    file = cellstr(strread(file,'%s'));
    
    %file = strsplit(file);

    %names = [path2 '/ipp*'];
    %name1 = dir(/mntdirect/_data_bm05_inhouse/seb/140312_colloids+Zeiss/Scans_0.45um_objects/Flower_Petalpath2);
    file1 = imread(file{1});
    file1 = file1(10:end,:);
end;
if myleefilter > 0, file1 = myerosion(file1, myleefilter); end; 
disp('  DONE');    fprintf('\n\n');

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
    [files] = open_seq(path1,ROI1);
    %-------------flat and darkfield image corection------------%
    flatim1 = []; darkim = [];
    if ~isempty(darkfieldImagePath) && exist(darkfieldImagePath,'dir')
        [darkheader, darkim] = pmedf_read(fullfile(darkfieldImagePath,darkfieldImageName));
        darkim = darkim(ROI1(1):ROI1(2),ROI1(3):ROI1(4));
        if rotim ~= 0, darkim = rot90(darkim,rotim); end;
    end;
    if ~isempty(FlatfieldImagePath) && exist(FlatfieldImagePath,'dir')
        [flatheader,flatim1] = pmedf_read(fullfile(FlatfieldImagePath,FlatfieldImageName1));
        flatim1 = flatim1(ROI1(1):ROI1(2),ROI1(3):ROI1(4));
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
        darkim = darkim(ROI1(1):ROI1(2),ROI1(3):ROI1(4));
        if rotim ~= 0, darkim = rot90(darkim,rotim); end;
    end;
    if ~isempty(FlatfieldImagePath) && exist(FlatfieldImagePath,'dir')
        flatim1 = imread(fullfile(FlatfieldImagePath,FlatfieldImageName1),'tif');
        flatim1 = flatim1(ROI1(1):ROI1(2),ROI1(3):ROI1(4));
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
    %files(k).data(ROI1(1):ROI1(2),ROI1(3):ROI1(4));%already cropped
    if myleefilter > 0, files(k).data = myerosion(files(k).data, myleefilter,'silent'); end; 
    if filtersz > 0 ,files(k).data = uint16(imfilter(files(k).data,ones(filtersz)./(filtersz).^2,'replicate','same'));end;
    if exist('rotim','var')   && (rotim ~= 0), files(k).data = rot90(files(k).data,rotim );end;
    if ~isempty(darkim), files(k).data = single(files(k).data) - darkim;end;
    if ~isempty(flatim1), files(k).data = single(files(k).data) ./ single(flatim1);end;
    if undersampl > 0,files(k).data = files(k).data(1:undersampl:end,1:undersampl:end); end;
    progmeter(k/nImages);
end;
if undersampl > 0, pixsize = pixsize *undersampl;end;

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
if strcmp(synchrotron , 'ESRF'),    [files] = open_seq(path2,ROI1);
elseif strcmp(synchrotron , 'DLS'),
   [files] = open_seq_dia(path2,ROI1);
   %-------------flat and darkfield image corection------------%
   if ~isempty(FlatfieldImagePath) && exist(FlatfieldImagePath,'dir')
        flatim2 = imread(fullfile(FlatfieldImagePath,FlatfieldImageName2),'tif');
        flatim2 = flatim2(ROI1(1):ROI1(2),ROI1(3):ROI1(4));
        if rotim > 0, flatim2 = rot90(flatim2,rotim); end;
        if myleefilter > 0, flatim2 = myerosion(flatim2, myleefilter); end;       
        if ~isempty(darkim), flatim2 = flatim2- darkim;end; 
   end;
end;
disp(' Flatfield images loaded');

% error checking
if ~isstruct(files) || isempty(files(1)), error('No FlatField Spictures in memory: path probably not correct');end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% reduce to ROI and process
for k = 1 : 1 : nImages, 
    files(k).data = files(k).data;
    %files(k).data(ROI1(1):ROI1(2),ROI1(3):ROI1(4));%already cropped
    if myleefilter > 0, files(k).data = myerosion(files(k).data, myleefilter,'silent'); end; 
    if filtersz > 0 ,files(k).data = uint16(imfilter(files(k).data,ones(filtersz)./(filtersz).^2,'replicate','same'));end;
    if exist('rotim','var')   && (rotim ~= 0), files(k).data = rot90(files(k).data,rotim );end;
    if ~isempty(darkim), files(k).data = single(files(k).data) - darkim;end;
    if ~isempty(flatim1), files(k).data = single(files(k).data) ./ single(flatim1);end;
    if undersampl > 0,files(k).data = files(k).data(1:undersampl:end,1:undersampl:end); end;
    progmeter(k/nImages);
end;
if undersampl > 0, pixsize = pixsize *undersampl;end;
figure(6)
subplot(1,2,2)
imagesc(files(1).data)
colorbar
colormap(gray)
title('First image ref')


% %% %%%%%%%%%%%%%%%%%%%%%%%%%%%% correct membrane position error%%%%%%%%%%%%55555
% if correctMembPo,
%         for k = 1:
%             
% end;



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
        imgr = files(pp).data;
        piece_stack_ref{k}(:,:,pp) = imgr(zone,:);
        piece_stack_sample{k} = uint16(piece_stack_sample{k});
        piece_stack_ref{k} = uint16(piece_stack_ref{k});
    end;
end;


stsample = cat(1,piece_stack_sample{:}); 
stref = cat(1,piece_stack_ref{:});
%clear piece_stack_sample piece_stack_ref
[m,n,r] = size(stsample);

return;