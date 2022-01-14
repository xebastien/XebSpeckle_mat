 %% calculate the displacment between two set of images from raster scans
% first version used for the paper 'X-ray multimodal imaging with..."
% this code is slow and an other one under dev at08/2013
function I = pixdelaycorr_findI0(stack_sample,stack_ref)
%edgeSizeVH,
%edgeSize.V = edgeSizeVH(1); edgeSize.H = edgeSizeVH(2);

if (sqrt(size(stack_sample,3))-round(sqrt(size(stack_sample,3))))>0,error('Non consistent number of images');end;
if size(stack_sample,3) ~= size(stack_ref,3), error('stacks dont have the same sie');end;

s = size(stack_sample);
nx = s(1);
ny = s(2);
nxy = nx*ny;
n = s(3);
nr = sqrt(n);% lateral size of the scan
h = fspecial('gaussian', 3);% 0 for no filtering

newssample = zeros(nr,nr,nxy);
newsref = zeros(nr,nr,nxy);

% images normalization
stack_sample = single(stack_sample);stack_ref = single(stack_ref);
% for pp = 1:size(stack_sample,3)
%     stack_sample(:,:,pp)  =  (stack_sample(:,:,pp) - mean(reshape(stack_sample(:,:,pp),1,[])))./std(reshape(stack_sample(:,:,pp),1,[]));
%     stack_ref(:,:,pp)     =  (stack_ref(:,:,pp)  -   mean(reshape(stack_ref(:,:,pp),1,[])))   ./std(reshape(stack_ref(:,:,pp),1,[])   );
% end;



for pp = 1:nr
    for qq = 1: nr
            newssample(qq,pp,:) = reshape(stack_sample(:,:,(pp-1).*nr+qq),[],1); 
            newsref(qq,pp,:) = reshape(stack_ref(:,:,(pp-1).*nr+qq),[],1);
    end;
end;
    
%thephase = zeros(2,nxy);
% edgeSize.V = 3;%max([edgeSize.V1 edgeSize.V2]);
% edgeSize.H = 3;%max([edgeSize.H1 edgeSize.H2]);
%Iav = [0 0];
%filtersz = 3;
samplingRate = 3;
 I = zeros(length(1:samplingRate:nxy),2);
% %figure(12)
% 
% 
%# crosscorrelation
edgeS = -1;

for pp = 1:samplingRate:nxy
        ns = newssample(:,:,pp);
        nf = newsref(:,:,pp);
        %ns = medfilt2(ns,[3 3]);nf = medfilt2(nf,[3 3]);
        ns = imfilter(ns,ones(3),'same','replicate');nf = imfilter(nf,ones(3),'same','replicate');
        ns([1:edgeS+1 (nr-edgeS):nr],:)=0;ns(:,[1:edgeS+1 (nr-edgeS):nr])=0;
        
        cor0 = conj(fft2(newssample(:,:,pp))).*fft2(nf);  
        cor3 = fftshift(ifft2( cor0));
        %cor2(floor(nr/2)+1,floor(nr/2)+1) = mean(cor2(:));
        %cor2(floor(nr/2):floor(nr/2)+1,floor(nr/2):floor(nr/2)+1) = mean(cor2(:));
        [~,I(floor(pp/samplingRate)+1,:),I(floor(pp/samplingRate)+1,2)] = max2(cor3);
%       [cor2,pic3] = normxcorr2_mexsub(newssample(6:end-5,6:end-5,pp),nf,'same');
%         
%         
%         figure(2)
% %         subplot(2,2,1)
% %         imagesc(((cor2)));
%         subplot(2,2,3)
%         imagesc((abs(ns)));
%         subplot(2,2,2)
%         imagesc((abs(nf)));
%         subplot(2,2,4)
%         imagesc(((cor3)));
        
end;
 Iav = I - repmat(floor(size(cor3)./2)+1,length(1:samplingRate:nxy),1);
I = median(Iav,1);
%disp(num2str(I));
%I = -[5 2];
%disp(num2str(I));
% figure(1)
% plot(Iav)
% drawnow


return;
