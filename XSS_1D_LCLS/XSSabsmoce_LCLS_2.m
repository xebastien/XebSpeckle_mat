%% input parameters
% hard parameters that you shouldnt touch unless you know why
edgeSizeStep = 45;   % larger  number of step disaplcment
edgeSizePix = 25 ;    % larger number of pixel displacment for the disp vector
bandwidth = 3;%must be smaller than edgeSizePix

%% ===================================================================================
%                                   find rigid displacement
% ====================================================================================
szspot = 30;
winvect = 3;        % width of the vector in 2D *2+1

winsize = 11;       % who far to probe and look for the correlation peak

[X,Y] = meshgrid(1:1:n1,1:1:m1);
file1 = single(files1(1).data);
ROIc = [1 size(file1,1) 1 size(file1,2)];
CentMassX = round((ROIc(3) + ROIc( 4))./2);
CentMassY = round((ROIc(1) + ROIc(2))./2);

avg1 = ones(size(files1(1).data));
avg2 = ones(size(files2(1).data));

for k = 1:nImages,
    files1(k).data(isnan(files1(k).data)) = 0;
    avg1 = avg1 + double(files1(k).data)./nImages;
    avg2 = avg2 + double(files2(k).data)./nImages;
end;
spot1n = double(files1(1).data)./avg1;        spot2n = double(files2(1).data)./avg2;

spot2n(isnan(spot2n)) = 0;
spot1n = spot1n(CentMassY-szspot:CentMassY+szspot,CentMassX-szspot:CentMassX+szspot);
% find the rigid translation0
[cor,A] = normxcorr2_mexsub(double(spot1n),spot2n,'valid');

[~,mcent(2)] = max(max(cor,[],1),[],2);
[~,mcent(1)] = max(max(cor,[],2),[],1);



szfc.offsetV = mcent(1) - CentMassY+szspot;
szfc.offsetH = mcent(2) - CentMassX+szspot;
disp(['Rigid translation of the beam :' num2str([szfc.offsetV szfc.offsetH]) ' V - H   pixels'])
return;
files2(1).data = circshift(files2(1).data,round(mcent(1)-size(files2(1).data,1)./2),1);
files2(1).data = circshift(files2(1).data,round(mcent(2)-size(files2(1).data,2)./2),2);



%% create new sub stack for fast processing with multicores
[m1,n1] = size(files1(1).data);
[m2,n2] = size(files2(1).data);

zonebounds1 = round(linspace(edgeSizePix,m1-edgeSizePix,subdiv+1));
zonebounds2 = round(linspace(edgeSizePix,m2-edgeSizePix,subdiv+1));
clear piece_stack_sample2 piece_stack_ref2
piece_stack_sample2 = cell(1,subdiv);
piece_stack_ref2 = cell(1,subdiv);
for k = 1 : subdiv
    zone1 = -edgeSizePix + zonebounds1(k) + 1:1:zonebounds1(k+1)+ edgeSizePix;
    zone2 = -edgeSizePix + zonebounds2(k) + 1:1:zonebounds2(k+1)+ edgeSizePix;
    piece_stack_sample2{k} = zeros(length(zone1) , n1, nImages );
    piece_stack_ref2{k} = zeros(length(zone2) , n2, nImages );
    for pp = 1 : nImages, 
        piece_stack_ref2{k}(:,:,pp) = files2(pp).data(zone2,:); 
        piece_stack_sample2{k}(:,:,pp) = files1(pp).data(zone1,:);
    end;
end;
%% ---------------------do the job  ---------------------------
elphase = cell(1,subdiv);
for k = 1: subdiv    
    elphase{k} = galPhase1D_onMEX(piece_stack_sample2{k},piece_stack_ref2{k},edgeSizeStep,edgeSizePix,bandwidth);
end;


elphase = cat(1,elphase{:});
