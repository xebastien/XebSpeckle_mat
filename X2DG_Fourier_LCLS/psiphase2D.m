% 2D shearing interferometry using 2D phase stepping
%
% by S.berujon 
% oct 2011
% from a stack of N^2 picture from a raster scan of a 2D interferometer,
% returns the 2 different phase map, the absorption map, and the 4
% scattering maps
%
%INPUT : 
% stack :size must be  m x n x g. it contains the 16 images from the scan
%
% OUTPUT : ampMap : the absorption image
% phaseMaps: is of size m x n x 2 . it contains the 2 phase maps (resp vert and horiz)
% scatMaps : is of size m x n 4. it contains the 4 darkfield images )
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ampMap,phaseMaps,scatMaps ] = psiphase2D( stack,freq,PadFact )  

%if ~((sqrt(size(stack,3))== 4) || (sqrt(size(stack,3)) == 8) || (sqrt(size(stack,3)) == 16)), error('stack has not the correct expected size(16)');end;


[nx,ny,n] = size(stack);
nxy = nx*ny;

nr = sqrt(n);
 
% Build pixel stack


newstack = zeros(nr,nr,nxy);
for pp = 1:1:nr
    for qq = 1: 1:nr
        newstack(qq,pp,:) =reshape(stack(:,:,(pp-1).*nr+qq),[],1); 
    end;
end;
    
% windowing
bell1 = repmat(hamming(nr),[1 nr]).* repmat(hamming(nr),[1 nr])';

for k = 1:nxy
    newstack(:,:,k) = (myerosion(newstack(:,:,k),3,3,'silent') - mean(reshape(newstack(:,:,k),1,[]))).*bell1;
end; %#ok<*NOSEL>
%     
    
% for kk = 1:200
%     figure(111)
%     imagesc(newstack(:,:,kk.*50))
%     colormap(gray)
%     drawnow
%     pause(.2)
% end;
% newstack2 = zeros(64,64, nxy);
% newstack2(1:nr,1:nr,:) = newstack;


newstack2 = zeros(nr*PadFact,nr*PadFact,nxy);
newstack2(1:nr,1:nr,:) = newstack;


ft = fft2(newstack2);

phaseMapH = angle(ft(round(PadFact*freq)+1,round(PadFact*freq)+1,:));
phaseMapV = angle(ft(round(PadFact*freq)+1,end-round(PadFact*freq),:));

phaseMapsV = squeeze(reshape(phaseMapV,[nx ny 1]));
phaseMapsH = squeeze(reshape(phaseMapH,[nx ny 1]));

phaseMaps = zeros(nx,ny,2);
phaseMaps(:,:,1) = phaseMapsV;
phaseMaps(:,:,2) = phaseMapsH;

ampMap = reshape(abs(ft(1,1,:)),[nx ny 1]);



visMapX = abs(  ft(1,round(PadFact*freq)+1,:)   )./abs(  ft(1,1,:)   );
visMapY = abs(  ft(round(PadFact*freq)+1,1,:)   )./abs(  ft(1,1,:)   );

visMapDX = abs(  ft(4,4,:)   )./abs(  ft(1,1,:)   );
visMapDY = abs(  ft(round(PadFact*freq)+1,round(PadFact*freq)+1,:)   )./abs(  ft(1,1,:)   );

scatMaps = zeros(nx,ny,4);
scatMaps(:,:,1) = squeeze(reshape(visMapX,[nx ny 1]));
scatMaps(:,:,2) = squeeze(reshape(visMapY,[nx ny 1]));
scatMaps(:,:,3) = squeeze(reshape(visMapDX,[nx ny 1]));
scatMaps(:,:,4) = squeeze(reshape(visMapDY,[nx ny 1]));