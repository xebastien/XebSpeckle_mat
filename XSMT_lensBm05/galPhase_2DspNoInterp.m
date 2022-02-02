 % calculate the displacment between two set of images from raster scans
% first version used for the paper 'X-ray multimodal imaging with..."
% this code is slow and an other one under dev at08/2013
function laphase = galPhase_2DspNoInterp(stack_sample,stack_ref,undersamp)

if (sqrt(size(stack_sample,3))-round(sqrt(size(stack_sample,3)))) > 0,
    error('Non consistent number of images');
end;

if (sqrt(size(stack_ref,3))-round(sqrt(size(stack_ref,3)))) > 0,
    error('Non consistent number of images');
end;

s = size(stack_sample);
r = size(stack_ref);

nx = s(1);
ny = s(2);
nxy = nx*ny;
ns = s(3);      nssr = sqrt(ns);
nr = r(3);      nrsr = sqrt(nr);


newssample = zeros(nssr,nssr,nxy);
newsref = zeros(nrsr,nrsr,nxy);

for pp = 1:nssr
    for qq = 1: nssr
        newssample(qq,pp,:) = reshape(stack_sample(:,:,(pp-1).*nssr+qq),[],1); 
    end;
end;
    
for pp = 1:nrsr
    for qq = 1: nrsr
        newsref(qq,pp,:) = reshape(stack_ref(:,:,(pp-1).*nrsr+qq),[],1); 
    end;
end;
% newssample(1:12,1:12,:);
% newsref(1:12,1:12,:);

thephase = zeros(2,nxy);
%absorp = zeros(1,nxy);
edgeSize = (nrsr - ((nssr-1)*undersamp+1))/2+1;
sizesmall = 1 + undersamp*(nssr-1);
for pp = 1:1:nxy


    plane3t = double(squeeze(newsref(:,:,pp)));
    plane2t = zeros(sizesmall,sizesmall,'double');
    plane2t(1:undersamp:end,1:undersamp:end) = double(squeeze(newssample(:,:,pp)));
    
    [cor,pic3] = normxcorr2_mexsub(plane2t, plane3t,'valid');
    
    
    plane2tr = zeros(sizesmall,sizesmall,'double');
    plane2tr(1:undersamp:end,1:undersamp:end) = plane3t(edgeSize:undersamp:end-edgeSize+1,edgeSize:undersamp:end-edgeSize+1);
    [cor2,pic4] = normxcorr2_mexsub(plane2tr, plane3t,'valid');
    
    
    thephase(1,pp) = pic3(4)-pic4(4);
    thephase(2,pp) = pic3(3)-pic4(3);
    
    if abs(pic3(4)) > edgeSize || abs(pic3(3)) > edgeSize, 
        [thephase(1,pp), thephase(2,pp), ~] = findpeak(cor,true );
        thephase(2,pp) = thephase(2,pp) - edgeSize;
        thephase(1,pp) = thephase(1,pp) - edgeSize;
    end;
    
    if isnan(pic3(4)) || isnan(pic3(3)),
        [thephase(1,pp), thephase(2,pp), ~] = findpeak(cor,true );
        thephase(2,pp) = thephase(2,pp) - edgeSize;
        thephase(1,pp) = thephase(1,pp) - edgeSize;
    end;
%     figure(1)
%     subplot(1,3,1)
%     imagesc(plane1t)
%     subplot(1,3,2)
%     imagesc(squeeze(newsref(:,:,pp)))
%     subplot(1,3,3)
%     imagesc(cor)

end;

laphase(:,:,1) = reshape(thephase(1,:),nx,ny);
laphase(:,:,2) = reshape(thephase(2,:),nx,ny);


return;
