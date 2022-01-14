% XSS2D with grating structure 
% propagation distance plan

%
% LOAD only the images
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% global varia-ble setting

% noPlate_checkboard = 1 ; Plate_checkboard = 2;  hexa_noPlate = 3;
% hexa_Plate = 4; pi2.85_noPlate = 5;
scanChoice = 2;

switch scanChoice
    case 1
        scanName = 'run_checkb4um_noPlate.h5';
        ROI1 = [480 1641 320 1470];
    case 2
        scanName = 'run_checkb4um_Plate.h5';
        ROI1 = [480 1641 320 1470];
    case 3
        scanName = 'run_hexa_NoPlate.h5';
        ROI1 = [480 1641 320 1470];
    case 4
        scanName = 'run_hexa_Plate.h5';  
        ROI1 = [480 1641 320+93 1470+100];
    case 5
        scanName = 'run_2.85pi_noPlate.h5';  
        ROI1 = [480 1641 320+93 1470+160];
    case 6
        scanName = 'run50_bigPitch14um.h5';  
        ROI1 = [480 1641 320+93 1470+160];
        
end;

synchrotron ='LCLS';% 'LCLS' or 'ESRF'or 'DLS'for the format of images

filePath  = 'E:\LCLS';          %
maxImages = 17490;
subdiv    = 19;    % chunck to use to avoid memory issues.

im_return   = 0;% return afet display of the first image
filtersz    = 0;% averaging filter size
myleefilter = 30;% for th1 erosion filter (morphological processing): 0 = none, value = filter threshold
undersampl  = 0; % positive for binning
interpsecim = 0;%(1.4338673051/1.20776665922)^(1);    %          %  pixel size 2 / (pixel size 1) 1.4184
rotim = -1;


% -------------------------dark field and flat field images--------------%
defaultPath = '/data/bm05/imaging/seb/xfel';% where to pick the first images


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    open the files of the sample
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fileName = fullfile(filePath,scanName);
    
fprintf('\n Starting loading pictures\n_____');
disp(fileName);
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
    clear det1
end;


% error checking
if ~isstruct(files) || isempty(files(1)), error('No 1pictures in memory: path probably not correct');end;
nImages = length(files);
if nImages > maxImages, nImages = maxImages;end;
disp(['Number of pictures of the sample = ' num2str(nImages) ' or ' num2str(sqrt(nImages)) ' ^2'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% reduce to ROI and process
for k = 1 : 1 : nImages, 
    if myleefilter > 0,       files(k).data = myerosion(files(k).data, myleefilter,'silent');         end; 
    if filtersz > 0 ,         files(k).data = uint16(imfilter(files(k).data,ones(filtersz)./(filtersz).^2,'replicate','same'));end;
    if exist('rotim','var'),  files(k).data = fliplr(rot90(files(k).data,rotim ));                    end;
    if ~isempty(darkim),      files(k).data = single(files(k).data) - darkim;                         end;
    if ~isempty(flatim1),     files(k).data = single(files(k).data) ./ single(flatim1);               end;
    if undersampl > 0,        files(k).data = files(k).data(1:undersampl:end,1:undersampl:end);       end;
    progmeter(k/nImages);
end;
%if undersampl > 0, pixsize = pixsize *undersampl;end;


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           reshape the data                           %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[m1,n1] = size(files(1).data);
disp('  DONE');

stsample = zeros(m1,n1,length(files));%*length(scanNumbers));

for k = 1 :length(files)
    stsample(:,:,k) = (single(files(k).data) - mean(single(files(k).data(:))))./std(single(files(k).data(:)));
end;
stsample = stsample(ROI1(1):ROI1(2),ROI1(3):ROI1(4),:);
%%%%  small plots of the first image
%clear piece_stack_sample piece_stack_ref
[m,n,nImages] = size(stsample);
figure(6)
subplot(1,2,1)
imagesc(files(1).data)
subplot(1,2,2)
imagesc(stsample(:,:,end))
colorbar
colormap(gray)
title('First image ref')
clear files m n m1 n1 
return;