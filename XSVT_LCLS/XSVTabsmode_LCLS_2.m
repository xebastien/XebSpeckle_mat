        %% ===================================================================================
%                                   find rigid displacement
% ====================================================================================
szspot = 30;    %to calculate the rigid tranlsation between the two images
winvect = 2;    % width of the vector in 2D 

winsize = 13;   % who far to probe and look for the correlation peak - larger is longer but better dynamic of the technique ( max grad calc)

ROI1 = [475 975 240 722  ].*2;
%ROI1 = [400 1550 580 1700];% for 1 trip one lens euXfel
[X,Y] = meshgrid(1:1:n1,1:1:m1);
% for the caculation of the mask later
file1 = sum(stsample(ROI1(1)+winsize+1:ROI1(2)-winsize-1,ROI1(3)+winsize+1:ROI1(4)-winsize-1,:),3);
ROIc = ROI1;
%% looking for the rigid tranlation between the two images. We correlate a single bis subset
CentMassX = round((ROIc(3) + ROIc( 4))./2);
CentMassY = round((ROIc(1) + ROIc(2))./2);

avg1 = ones(size(files1(1).data));
avg2 = ones(size(files2(1).data));

for k = 1:nImages
    files1(k).data(isnan(files1(k).data)) = 0;
    avg1 = avg1 + double(files1(k).data)./nImages;
    avg2 = avg2 + double(files2(k).data)./nImages;
end;
spot1n = double(files1(1).data)./avg1;        spot2n = double(files2(1).data)./avg2;

spot2n(isnan(spot2n)) = 0;
spot1n = spot1n(CentMassY-szspot:CentMassY+szspot,CentMassX-szspot:CentMassX+szspot);
% find the rigid translation0
[cor,A] = normxcorr2_mexsub(double(spot1n),spot2n,'valid');

%[~,mcent(1),mcent(2)] = max2(cor);
[~,mcent(2)] = max(max(cor,[],1),[],2);
[~,mcent(1)] = max(max(cor,[],2),[],1);


szfc.offsetV = mcent(1) - CentMassY + szspot;
szfc.offsetH = mcent(2) - CentMassX + szspot;
disp(['Rigid translation of the beam :' num2str([szfc.offsetV szfc.offsetH]) ' V - H   pixels'])
% if you see a nice little peak, you won. Otherwise you messed up somewhere
figure(1)
subplot(2,2,1)
imagesc(files1(1).data)
subplot(2,2,2)
imagesc(files2(1).data)
colormap jet
subplot(2,2,3)
imagesc(cor)
subplot(2,2,4)
imagesc(spot1n)

%% recrop the stack to be centered on the beam for each stack

stsample1 = stsample(ROI1(1):ROI1(2),ROI1(3):ROI1(4),:);
stref1 = stref(ROI1(1)+szfc.offsetV:ROI1(2)+szfc.offsetV,ROI1(3)+szfc.offsetH:ROI1(4)+szfc.offsetH,:);
stsample1 = single(stsample1 );
stref1 = single(stref1);

stref1    = stref1./repmat(mean(stref1,3),[1 1 nImages]);
stsample1 = stsample1./repmat(mean(stsample1,3),[1 1 nImages]);

%% we enlarge the 3D matrix along the third dimension by putting also the
% value of the neighboring pixels in the vectors
if winvect<1,winvect = 1;end;
strefBig = [];stsampBig = [];% these are the artifically enlarge staks we are going to use

for kv = 1:1:winvect
    for kh = 1:1:winvect
        strefBig = cat(3,strefBig,stref1(kv : end - winvect + kv,kh : end - winvect + kh,:)); 
        stsampBig = cat(3,stsampBig,stsample1(kv : end - winvect + kv,kh : end - winvect + kh,:));
    end;
end;

stref1 = strefBig;
stsample1 = stsampBig;

absorption = mean(stsampBig,3)./mean(strefBig,3);
absorption = absorption(winsize+1:end-winsize-1,winsize+1:end-winsize-1);
file1 = absorption;

clear stsampBig strefBig
%% calculate the vector dispalcment

stsample1 = (stsample1 - repmat(mean(stsample1,3),[1 1 size(stsample1,3)]))./repmat(std(stsample1,[],3),[1 1 size(stsample1,3)]);
stref1 = (stref1 - repmat(mean(stref1,3),[1 1 size(stref1,3)]))./repmat(std(stref1,[],3),[1 1 size(stref1,3)]);
stref1(isnan(stref1)) = 1;stsample1(isnan(stsample1)) = 1;
nCores = 20*3-5;


[m,n,p] = size(stref1);
zonebounds = round(linspace(1+winsize,m-winsize,nCores+2));
zonebounds = [zonebounds(1:end-2) zonebounds(end)]; 
out1 = cell(1,nCores);out2 = cell(1,nCores);out3 = cell(1,nCores);sts = cell(1,nCores);str = cell(1,nCores);
for k =nCores:-1:1
    sts{k} = double(stsample1(zonebounds(k)-winsize:zonebounds(k+1)+winsize,:,:));
    str{k} = double(stref1(zonebounds(k)-winsize:zonebounds(k+1)+winsize,:,:));
end;

parfor k = 1:1:nCores
    out1{k} = trackvect_mex(str{k},sts{k},int8(winsize));
%     out2{k} = trackvect_mex(sts{k},sts{k},int8(winsize));
%     out3{k} = trackvect_mex(str{k},str{k},int8(winsize));
    disp([num2str(k) 'Done']);
end

out1 = cat(1,out1{:});
% out3 = cat(1,out3{:});
% out2 = cat(1,out2{:});
disp('Vector displacement calculated');

% out4 = out1 - out2./2 - out3./2;