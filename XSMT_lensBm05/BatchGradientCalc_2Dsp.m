% batch processing for lens.
% calculate the wavefront gradient from XSS2S sparse mesh
%
%use PARTA AND PART B subfoolders
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                          
% if you want to use matlab2013a
% call for matlabpool to run with multiple cores (do it once at start):

% matlabpool(12)
% ------------------folders paths----------------%

rootFolder = '/data/bm05/inhouse/thomas/171211lenses/speckles/';
startName = 'lens' ;% how start the name of the folders you want to precess in this rootFolder
% name of the reference folders
refFolder1 = 'refs';
refFolder2 = 'refs';%'RXOPTICS_1D_ref_end';
partName = '*_2Dmesh_';
partNameRef1 = '*Start_2Dmesh_';
partNameRef2 = '*End_mesh_';
% ------------------export folder----------------%
% check that you have the writing permission in this folder
exportFolder = '/data/bm05/inhouse/thomas/171211lenses/processing/lens_2Dsp/';
%%%%%%%%%%                  filters                 %%%%%%%%%%%%%%%%%%%%%%%
rotim = -1;          % rot90  the pictures (Frelon at ESRF = -1 )
filtersz = 0;        % averaging filter size
myleefilter = 0;     % for the erosion filter (morphological processing): 0 = none, value = filter threshold (salt and pepper noise)
undersampl = 0;      % positive >1  for binning
h = 0;               % fspecial('gaussian', 3);    % 0 for no filtering
                     % miscelleanous--------------------------------------
maxImages = 5555;    % number of max images - ignore images beyond this number
subdiv = 20;         % number chuncks to create to avoid memory issues and also compute in parallel 

% hard parameters that you shouldnt touch unless you know why
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    Summary file                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



exportFileSummary = fullfile(exportFolder,'r1ameters_1step.txt');

dlmwrite(exportFileSummary, 'Processing parameters','delimiter', '' );
dlmwrite(exportFileSummary, ['rotim = ' num2str(rotim)], '-append','delimiter', '' );
dlmwrite(exportFileSummary, ['filtersz = ' num2str(filtersz)], '-append','delimiter', '' );
dlmwrite(exportFileSummary, ['myleefilter = ' num2str(myleefilter)], '-append','delimiter', '' );
dlmwrite(exportFileSummary, ['undersampl = ' num2str(undersampl)], '-append','delimiter', '' );
dlmwrite(exportFileSummary, 'special filter = ' , '-append','delimiter', '' );
dlmwrite(exportFileSummary, h(:)', '-append','delimiter', '\t'  ); %'
dlmwrite(exportFileSummary, ['subdiv = ' num2str(subdiv)], '-append','delimiter', '' );
dlmwrite(exportFileSummary, ['winsize = ' num2str(winsize)], '-append','delimiter', '' );
dlmwrite(exportFileSummary, 'Processed:', '-append','delimiter', '' );



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    select the folders to process                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subFolders = dir([rootFolder '/' startName '*']);

%subFolders = subFolders(3:end);
for k = length(subFolders):-1:1,
    if strcmp(subFolders(k).name,refFolder1) || strcmp(subFolders(k).name,refFolder2),
        subFolders(k) = [];
    end
end;

%subFolders = nestedSortStruct(subFolders, 'datenum');

for kname1 = 6:1:(length(subFolders)),
kname  = kname1;
% if       mod(kname1,2) && onlyPartA,  partName = 'partA';%partName = 'partA';
% elseif  ~mod(kname1,2) && onlyPartA,  continue;
% elseif   mod(kname1,2) && ~onlyPartA, partName = 'partA';
% else                                  partName = 'partB';end;
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    open the files of the sample                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\n Starting loading objects pictures from \n_____');
disp(fullfile(rootFolder,subFolders(kname).name,partName));fprintf('\n\n ');
[files] = open_seq(fullfile(rootFolder,subFolders(kname).name),'full',partName);

% error checking
if ~isstruct(files) || isempty(files(1)), error('No 1pictures in memory: path probably not correct');end;
nImages = length(files);
if nImages > maxImages,     nImages = maxImages;        end;
disp(['Number of pictures of the sample = ' num2str(nImages)])

%--------------------------------------------------------------------------
% preprocess
for k = 1 : 1 : nImages, 
    if myleefilter > 0, files(k).data = myerosion(files(k).data, myleefilter,'silent'); end; 
    %if exist('rotim','var')   && (rotim ~= 0), files(k).data = rot90(files(k).data,rotim );end;
    if filtersz > 0 ,files(k).data = uint16(imfilter(files(k).data,ones(filtersz)./(filtersz).^2,'replicate','same'));end;
    if undersampl > 0,imfiltered = imfilter(files(k).data,ones(3)./9,'same'); files(k).data= imfiltered(1:undersampl:end,1:undersampl:end); end;
    progmeter(k/nImages);
    if h(1) ~= 0, files(k).data = imfilter(files(k).data,h,'same');    end;
end;


[m1, n1] = size(files(1).data);
% _______________________Plot the first image______________________________
disp('***************');
[m, n] = size(files(1).data);
% figure(6)
% subplot(1,2,1)
% imagesc(rot90(files(1).data,rotim))
% colorbar
% colormap(gray)
% title('First object image')
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                  reference image loading                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

if (kname ) <= (length(subFolders)), refFolder = refFolder1;partNameRef = partNameRef1;
else refFolder = refFolder2;partNameRef = partNameRef2;
end;


disp(['Loading ref images from  ' refFolder]);fprintf('\n\n ');
filesSample = files;% save for the 
clear files
% open the files from the path
fprintf('\n Starting loading references pictures from o\n_____');
disp(refFolder);
[files] = open_seq(fullfile(rootFolder,refFolder),'full',partNameRef);

disp('References images loaded');

% error checking
if ~isstruct(files) || isempty(files(1)), error('No FlatField Spictures in memory: path probably not correct');end;
%---------------------------------------------------------------------------
%%% preprocess
for k = 1 : 1 : length(files), 
    if myleefilter > 0, files(k).data = myerosion(files(k).data, myleefilter,'silent'); end; 
    %if exist('rotim','var')   && (rotim ~= 0), files(k).data = rot90(files(k).data,rotim );end;
    if filtersz > 0 ,files(k).data = uint16(imfilter(files(k).data,ones(filtersz)./(filtersz).^2,'replicate','same'));end;
    if undersampl > 0,imfiltered = imfilter(files(k).data,ones(3)./9,'same'); files(k).data= imfiltered(1:undersampl:end,1:undersampl:end); end;
    progmeter(k/nImages);
    if h(1) ~= 0,files(k).data = imfilter(files(k).data,h,'same'); end;
end;

% figure(6)
% subplot(1,2,2)
% imagesc(files(1).data)
% colorbar
% colormap(gray)
% title('First image ref')

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              pad with zeros missing images       and normalization                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


nImagesR = length(files);

% stref = zeros(m,n,nImagesR);
% stsamp = zeros(m,n,nImages);
% for k= 1 :length(files),    stref(:,:,k) = single(files(k).data);   end;
% for k= 1 :length(filesSample),    stsamp(:,:,k) = single(filesSample(k).data);   end; 
%     
% absorp = mean( stsamp,3) ./ mean( stref,3);
% scattering = (std(stsamp,[],3)./std(stref,[],3))./absorp;
% 
% clear stsamp stref
avr = zeros(1,length(files));str = zeros(1,length(files));
for  k = 1:length(files), 
    %files(k).data = (single(files(k).data) - mean(single(files(k).data(:))) )./std(single(files(k).data(:)));
    avr(k) = mean(single(files(k).data(:)));
    str(k) = std(single(files(k).data(:)));
    %filesSample2(k).data = files(k).data.*0;
end;

avs = zeros(1,length(files));sts = zeros(1,length(files));
for  k = 1:length(filesSample), 
    %filesSample(k).data = (single(filesSample(k).data) - mean(single(filesSample(k).data(:))))./std(single(filesSample(k).data(:)));
    avs(k) = mean(single(filesSample(k).data(:)));
    sts(k) =  std(single(filesSample(k).data(:)));
end;


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            reshape the data                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

zonebounds = round(linspace(1,m1+1,subdiv+1));
piece_stack_sample = cell(1,subdiv);
piece_stack_ref = cell(1,subdiv);
for k = 1 :subdiv
    zone = zonebounds(k):1:zonebounds(k+1)-1;
    piece_stack_sample{k} = zeros(length(zone) , n1, nImages);
    piece_stack_ref{k} = zeros(length(zone) , n1, nImagesR );
    for pp = 1 : nImages, 
        imgs = filesSample(pp).data;
        piece_stack_sample{k}(:,:,pp) = imgs(zone,:);
    end;
    for pp = 1 : nImagesR, 
        imgr = files(pp).data;
        piece_stack_ref{k}(:,:,pp) = imgr(zone,:);
    end;
end;
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            track the vectors                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%do the magic
tic;

out = cell(1,subdiv);           outref = cell(1,subdiv);
absk = cell(1,subdiv);          scat = cell(1,subdiv);



parfor k =1:subdiv
    absk{k} = mean(double(piece_stack_sample{k}),3)./mean(double(piece_stack_ref{k}),3);
    scat{k} = std(double(piece_stack_sample{k}),[],3)./std(double(piece_stack_ref{k}),[],3)./absk{k};
    
    for pp = 1 : nImagesR
        piece_stack_ref{k}(:,:,pp)  = (piece_stack_ref{k}(:,:,pp) - avr(pp))./str(pp);
    end;
    for pp = 1 : nImages
        piece_stack_sample{k}(:,:,pp)  = (piece_stack_sample{k}(:,:,pp) - avs(pp))./sts(pp);
    end
    
    out{k} = galPhase_2DspNoInterp(piece_stack_sample{k},piece_stack_ref{k},3);
%     out2{k} = galPhase_2DspNoInterp(piece_stack_sample{k}(2:end-1,1:end-2,:),piece_stack_ref{k}(2:end-1,2:end-1,:),3);
%     out3{k} = galPhase_2DspNoInterp(piece_stack_sample{k}(1:end-2,2:end-1,:),piece_stack_ref{k}(2:end-1,2:end-1,:),3);
%     out4{k} = galPhase_2DspNoInterp(piece_stack_sample{k}(3:end  ,2:end-1,:)  ,piece_stack_ref{k}(2:end-1,2:end-1,:),3);
%     out5{k} = galPhase_2DspNoInterp(piece_stack_sample{k}(2:end-1,3:end,:)  ,piece_stack_ref{k}(2:end-1,2:end-1,:),3);

    
    % the little mex script that will calculate the displacement of each stsample vect relatively to the ref stacks vectors
    % take a while

end;
% out2 = cat(1,out2{:}) ;
% out3 = cat(1,out3{:}) ;
% out4 = cat(1,out4{:}) ;
% out5 = cat(1,out5{:}) ;

clear thephase
thephase(:,:,1:2) = cat(1,out{:}) ;
thephase(:,:,3) = cat(1,absk{:}) ;
thephase(:,:,4) = cat(1,scat{:}) ;
%- 1./2.*cat(2,outref{:}) - 1./2.*cat(2,outsamp{:});

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            export result matrix                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
exportFileName = fullfile(exportFolder,[subFolders(kname).name '_' regexprep(partName,'*','') '.mat']);
save(exportFileName,'thephase');

dlmwrite(exportFileSummary, exportFileName , '-append','delimiter', '' );
disp(['Done with ' subFolders(kname).name]);

end;

return;
