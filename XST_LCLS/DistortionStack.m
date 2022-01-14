% Calcualte a detector distortion from a stack of images taken following a 2d
% constant step mesh.

% S. Berujon 2017

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


winsize = 3;
folderImages = '/data/bm05/inhouse/Ruxandra/170220_eucall_3rd/detector_distortion_horz/';
ndiff = 1;


files = open_seq(folderImages);
nImages = length(files)/2;


for k = 1:nImages,
    files(k).data = rot90(files(k).data,-1);
    files(k) .data = (single(files(k).data) - mean(single(files(k).data(:))))./std(single(files(k).data(:)));
end;
[m1,n1] = size(files(1).data);


stackIm = zeros(m1,n1,nImages);

for k = 1:nImages,
    stackIm(:,:,k) = files(k).data;
end;

% end loading images
%%  Calculate rigid moitn

resol = 35;    szfc.sq_sz_largeV = 65;     szfc.sq_sz_largeH = 65;

spot1n = double(files(ndiff).data);spot2n = double(files(2).data);
spot1n(1:floor(m1/3),:) = 0;        spot1n(:,1:floor(n1/3)) = 0;
spot1n(end-floor(m1/3):end,:) = 0;  spot1n(:,end-floor(n1/3):end) = 0;

% find the rigid translation0
cor = fftshift(ifft2( conj(fft2(spot1n)).*fft2(spot2n,size(spot1n,1),size(spot1n,2)) ));
[~,mcent(1),mcent(2)] = max2(cor);
offV = mcent(1) - floor(size(cor,1)./2)-1;
offH = mcent(2) - floor(size(cor,2)./2)-1;
disp(['Rigid translation of the beam :' num2str([offV offH]) ' V - H resp.   pixels'])


stackIm = cat(3,stackIm(2:end-1,2:end-1,:),stackIm(1:end-2,2:end-1,:),stackIm(2:end-1,1:end-2,:),stackIm(3:end,2:end-1,:),stackIm(3:end,2:end-1,:));

%% calculate vector displacment
s1 = double(stackIm(offV+1+winsize:end-offV-winsize,offH+1+winsize:end-offH-winsize,1:end-ndiff));
s2 = double(stackIm(offV+1+winsize + offV:end-offV-winsize + offV,offH + 1+winsize + offH:end-offH-winsize + offH,ndiff+1:end));

[m2,n2] = size(s1);



s1 = (s1 - repmat(mean(s1,3),[1 1 (size(s1,3))]))./repmat(std(s1,[],3),[1 1 (size(s1,3))]);
s2 = (s2 - repmat(mean(s2,3),[1 1 (size(s1,3))]))./repmat(std(s2,[],3),[1 1 (size(s1,3))]);

out = trackvect_mex(s1,s2,int8(winsize));

return;


    