% XSS2D with grating structure 
% propagation distance plan

%
% LOAD only the images
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% global varia-ble setting

close all

synchrotron ='LCLS';% 'LCLS' or 'ESRF'or 'DLS'for the format of images

% Step is always 0.4 um

% Diamond checkboar pith d=2.8, d_checker = 4 um, pi/2@8keV h = 3.6 um
% Silicon d = 2.66, d_checker =3.76, type pi@9keV, h~67 um

% diamond hexagonal 1 pitch = 4 um, d_hexa = um, pi/2@8keV, h=3.6um
% diamond hexagonal 2 pitch = 16 um, d_hexa = um, pi/2@8keV, h=3.6um


% % not a square number of images -> 16x26
% fileNamePref = '/data/bm05/imaging/seb/LCLS/speckle/run';          %
% fileNameSuff = '_det2.h5';                                         % 
% scanNumbers = [37:52];

% Phase Plate +  Diamond Checkboard 4um
% pattern residual
fileNamePref = '/data/bm05/imaging/seb/LCLS/speckle/run';          %
fileNameSuff = '_det2.h5';                                         % 
scanNumbers = [69:92 94:95];

% NoPhase Plate +  Diamond Checkboard 4um
% pattern residual
fileNamePref = '/data/bm05/imaging/seb/LCLS/speckle/run';          %
fileNameSuff = '_det2.h5';                                         % 
scanNumbers = [97:106 108:110 112:123 124];
% 
% % NoPhase Plate + Hexagonal 4 um pi/2 phase grating
% fileNamePref = '/data/bm05/imaging/seb/LCLS/speckle/run';          %
% fileNameSuff = '_det2.h5';                                         % 
% scanNumbers = [126:151];
% % 
% % Phase Plate + Hexagonal 4 um pi/2 phase grating
% fileNamePref = '/data/bm05/imaging/seb/LCLS/speckle/run';          %
% fileNameSuff = '_det2.h5';                                         % 
% scanNumbers = [156: 165 167 169:183];

% % NoPhase Plate + 2.85 um pi phase grating
% fileNamePref = '/data/bm05/imaging/seb/LCLS/speckle/run';          %
% fileNameSuff = '_det2.h5';                                         % 
% scanNumbers = [205:208 210:231];
%%%%%%%%%%%%%%%%%%%
maxImages = 17490;
subdiv = 19;    % chunck to use to avoid memory issues.

im_return   = 0;% return afet display of the first image
filtersz    = 0;% averaging filter size
myleefilter = 50;% for th1 erosion filter (morphological processing): 0 = none, value = filter threshold
undersampl  = 0; % positive for binning
interpsecim = 0;%(1.4338673051/1.20776665922)^(1);    %          %  pixel size 2 / (pixel size 1) 1.4184
rotim = -1;


% -------------------------dark field and flat field images--------------%
defaultPath = '/data/bm05/imaging/seb/xfel';% where to pick the first images


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    open the files of the sample
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
files = cell(1,length(scanNumbers));


for kScans = 1:length(scanNumbers)
    fileName = [fileNamePref num2str(scanNumbers(kScans)) fileNameSuff];
    
    fprintf('\n Starting loading pictures\n_____');
    disp(fileName);
    flatim1 = []; darkim = [];
    if strcmp(synchrotron , 'ESRF')    
        [files{kScans}] = open_seq(path1);%
    elseif strcmp(synchrotron , 'LCLS')    
        det1 = h5read(fileName,'/avg');
        for k = 1:size(det1,3),
            files{kScans}(k).data = det1(:,:,k);
            files{kScans}(k).data(isnan(files{kScans}(k).data)) = 0;
            files{kScans}(k).data(files{kScans}(k).data < 0) = 0;
        end;      
    end;


    % error checking
    if ~isstruct(files{kScans}) || isempty(files{kScans}(1)), error('No 1pictures in memory: path probably not correct');end;
    nImages = length(files{kScans});
    if nImages > maxImages, nImages = maxImages;end;
    disp(['Number of pictures of the sample = ' num2str(nImages) ' or ' num2str(sqrt(nImages)) ' ^2'])

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% reduce to ROI and process
    for k = 1 : 1 : nImages, 
        files{kScans}(k).data = files{kScans}(k).data;
        if myleefilter > 0, files{kScans}(k).data = myerosion(files{kScans}(k).data, myleefilter,'silent'); end; 
        if filtersz > 0 ,files{kScans}(k).data = uint16(imfilter(files{kScans}(k).data,ones(filtersz)./(filtersz).^2,'replicate','same'));end;
        if exist('rotim','var')   && (rotim ~= 0), files{kScans}(k).data = fliplr(rot90(files{kScans}(k).data,rotim ));end;
        if ~isempty(darkim), files{kScans}(k).data = single(files{kScans}(k).data) - darkim;end;
        if ~isempty(flatim1), files{kScans}(k).data = single(files{kScans}(k).data) ./ single(flatim1);end;
        if undersampl > 0,files{kScans}(k).data = files{kScans}(k).data(1:undersampl:end,1:undersampl:end); end;
        progmeter(k/nImages);
    end;
    %if undersampl > 0, pixsize = pixsize *undersampl;end;


    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                           reshape the data                           %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [m1,n1] = size(files{kScans}(1).data);
        disp('  DONE');
end;



stsample = zeros(m1,n1,length(files{1}).^2);%*length(scanNumbers));

for kScans = 1:length(scanNumbers)
    d3index = (kScans-1)*length(files{1})+1;
    for k = 1 :length(files{kScans})
        stsample(:,:,d3index) = (single(files{kScans}(k).data) - mean(single(files{kScans}(k).data(:))))./std(single(files{kScans}(k).data(:)));
        d3index = d3index +1;
    end;
end;

%%%%  small plots of the first image
%clear piece_stack_sample piece_stack_ref
[m,n,nImages] = size(stsample);
figure(6)
subplot(1,2,1)
imagesc(files{1}(1).data)
subplot(1,2,2)
imagesc(files{end}(end).data)
colorbar
colormap(gray)
title('First image ref')
return;