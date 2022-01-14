% XSVT in absolute mode from two stacks of images taken at two different
% propagation distance plan

%
% LOAD only the images
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% global varia-ble setting

close all

synchrotron ='LCLS';% 'LCLS' or 'ESRF'or 'DLS'for the format of images
fileName = '/data/bm05/imaging/seb/LCLS/speckle/run236_det0.h5';        % with phase plate
fileName = '/data/bm05/imaging/seb/LCLS/speckle/run13_det0.h5';         % no phase plate
%%%%%%%%%%%%%%%%%%%
maxImages = 17490;
subdiv = 19;%144*2;% chunck to use to avoid memory issues.

filtersz = 0;       % averaging filter size
myleefilter = 10;   % for th1 erosion filter (morphological processing): 0 = none, value = filter threshold
undersampl = 0; % positive for binning
interpsecim = 2.62/2;%(1.4338673051/1.20776665922)^(1);    %          %  pixel size 2 / (pixel size 1) 1.4184
rotim = -1;
% LCLS pixel size pix1 = 1.44276, pix2 = 1.44978;
% EuXFEL pixel size pi1 = 1.4338673051 pix2 = 1.20776665922

% -------------------------dark field and flat field images--------------%
defaultPath = '/data/bm05/imaging/seb/xfel';% where to pick the first images

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    open the files of the sample
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\n Starting loading pictures\n_____');
flatim1 = []; darkim = [];
if strcmp(synchrotron , 'ESRF')    
    [files] = open_seq(path1);%
elseif strcmp(synchrotron , 'LCLS')    
    det1 = h5read(fileName,'/avg');
    for k = 1:size(det1,3),
        files(k).data = det1(:,:,k);
        files(k).data(isnan(files(k).data)) = 0;
        files(k).data(files(k).data < 0) = 0;
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
    if myleefilter > 0, files(k).data = myerosion(files(k).data, myleefilter,'silent'); end; 
    if filtersz > 0 ,files(k).data = uint16(imfilter(files(k).data,ones(filtersz)./(filtersz).^2,'replicate','same'));end;
    if exist('rotim','var')   && (rotim ~= 0), files(k).data = fliplr(rot90(files(k).data,rotim ));end;
    if ~isempty(darkim), files(k).data = single(files(k).data) - darkim;end;
    if ~isempty(flatim1), files(k).data = single(files(k).data) ./ single(flatim1);end;
    if undersampl > 0,files(k).data = files(k).data(1:undersampl:end,1:undersampl:end); end;
    progmeter(k/nImages);
end;
if undersampl > 0, pixsize = pixsize *undersampl;end;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                  det2 image loading
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

disp(fileName);
files1 = files;clear files
% open the files from the path
fprintf('\n Starting loading Flatfield pictures\n_____');
if strcmp(synchrotron , 'ESRF'),    
    [files] = open_seq(path2);
elseif strcmp(synchrotron , 'LCLS')    
    det2 = h5read(fileName,'/avg2');
    for k = 1:size(det1,3),
        files(k).data = det2(:,:,k);
        files(k).data(isnan(files(k).data)) = 0;
        files(k).data(files(k).data < 0) = 0;
    end;      
   
end;
disp(' Flatfield images loaded');
% error checking
if ~isstruct(files) || isempty(files(1)), error('No FlatField Spictures in memory: path probably not correct');end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[m1 ,n1] = size(files1(1).data); 
[XI,YI] = meshgrid(linspace(1,n1,round(n1/interpsecim)),linspace(1,m1,round(m1/interpsecim)));

%%% reduce to ROI and process
for k = 1 : 1 : nImages, 
    files(k).data = files(k).data;
    if myleefilter > 0, files(k).data = myerosion(files(k).data, myleefilter,'silent'); end; 
    if filtersz > 0 ,files(k).data = uint16(imfilter(files(k).data,ones(filtersz)./(filtersz).^2,'replicate','same'));end;
    if exist('rotim','var')   && (rotim ~= 0), files(k).data = fliplr(rot90(files(k).data,rotim ));end;
    if ~isempty(darkim), files(k).data = single(files(k).data) - darkim;end;
    if ~isempty(flatim1), files(k).data = single(files(k).data) ./ single(flatim1);end;
    if undersampl > 0,files(k).data = files(k).data(1:undersampl:end,1:undersampl:end); end;
    if interpsecim > 0, files(k).data = interp2(single(files(k).data),XI,YI,'linear'); end
    %files(k).data = imrotate(single(files(k).data),-3.5,'bilinear');
    progmeter(k/nImages);
end
if undersampl > 0, pixsize = pixsize *undersampl;end
files2 = files;clear files


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           reshape the data                           %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[m1,n1] = size(files1(1).data);
[m2,n2] = size(files2(1).data);


stsample = zeros(m1,n1,nImages);        stref = zeros(m2,n2,nImages);
for k = 1 :nImages
    stsample(:,:,k) = single(files1(k).data);
    stref(:,:,k)    = single(files2(k).data);
end
avsample = mean(stsample,3);    avref = mean(stref,3);
stsample = stsample./repmat(avsample,[1 1 nImages]);
stref    = stref   ./repmat(avref   ,[1 1 nImages]);
%%%%  small plots of the first image
%clear piece_stack_sample piece_stack_ref
[m,n,r] = size(stsample);
figure(6)
subplot(1,2,1)
imagesc(files1(1).data)
subplot(1,2,2)
imagesc(files2(1).data)
colorbar
colormap(gray)
title('First image ref')
return;