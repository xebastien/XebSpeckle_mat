 %% calculate the displacment between two set of images from raster scans
% first version used for the paper 'X-ray multimodal imaging with..."
% this code is slow and an other one under dev at08/2013
function [delay] = pixdelaycorr2(stack_sample,stack_ref,I)
%edgeSizeVH,
%edgeSize.V = edgeSizeVH(1); edgeSize.H = edgeSizeVH(2);

    minva = 5;
    minvb = minva;
    
if I(1) == 0, minva = 2;end ;  
if I(2) == 0, minvb = 2;end  ;    
    
if (sqrt(size(stack_sample,3))-round(sqrt(size(stack_sample,3))))>0,error('Non consistent number of images');end;
if size(stack_sample,3) ~= size(stack_ref,3), error('stacks dont have the same sie');end;

s = size(stack_sample);
nx = s(1);
ny = s(2);
nxy = nx*ny;
n = s(3);
nr = sqrt(n);% lateral size of the scan
h = 0;%fspecial('gaussian', 3);% 0 for no filtering

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
    
thephase = zeros(2,nxy);

for pp = 1:1:nxy
    
    
    if h(1) ~= 0,
%         newssample(:,:,pp) = medfilt2(newssample(:,:,pp),[3 3]);
%         newsref(:,:,pp) = medfilt2(newsref(:,:,pp),[3 3]);
        
        newssample(:,:,pp) = imfilter(newssample(:,:,pp),h,'same'); 
        newsref(:,:,pp) = imfilter(newsref(:,:,pp),h,'same'); 
    end;

    a = minva;              c = minvb; 
    b = nr - minva+1;       d = nr - minvb +1;
    
    a = a  - I(1);          b = b  - I(1);
    c = c  - I(2);          d = d  - I(2);
    
    if a < 1, a = 1; end;
    if c < 1, c = 1; end;
    if b > nr, b = nr; end;
    if d > nr, d = nr; end;
    
    plane1t = squeeze(newssample(a:b,c:d,pp));
    %plane1t([1:a b:nr],:) = 0;plane1t(:,[1:c d:nr]) = 0;
    
    
    if ~isempty(plane1t)
    [cor,pic3] = normxcorr2_mexsub(plane1t,squeeze(newsref(:,:,pp)),'same');

    thephase(1,pp) = pic3(4)-1 -c+1*minvb;%+ (c-minva)/2;
    thephase(2,pp) = pic3(3)-1 -a+1*minva;%+ (a-minva) /2;

    
    if isnan(pic3(4)) || isnan(pic3(3)),
        [thephase(1,pp), thephase(2,pp), ~] = findpeak(cor,true );
        thephase(2,pp) = thephase(2,pp) - floor(size(cor,1)/2)-1  - c + 2*minvb;
        thephase(1,pp) = thephase(1,pp) - floor(size(cor,2)/2)-1  - a + 2*minva;
    end;
%     figure(1)
%     subplot(2,2,1)
%     imagesc(squeeze(newssample(:,:,pp)))
%     subplot(2,2,2)
%     imagesc(squeeze(newsref(:,:,pp)))
%     subplot(2,2,3)
%     imagesc(cor)
%     imagesc(plane1t)
%     subplot(2,2,4)
%     imagesc(cor)
%     title(num2str(I))
%     disp(num2str(pic3'))%fliplr(round(thephase(:,pp))') - I
    end;
end;

delay(:,:,1) = reshape(thephase(1,:),nx,ny);
delay(:,:,2) = reshape(thephase(2,:),nx,ny);
return;
