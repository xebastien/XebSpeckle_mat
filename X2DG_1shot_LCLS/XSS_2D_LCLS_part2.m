% ===============================================================
pixreach.V  = 3; % number of interpixel distance to correlate
pixreach.H  = 3;
scanDimOrd  = 1; % )How was made the mesh scan: 1st dim is vert or horizontal ? boolean 1 or 0
PitchG      = 4;     StepS = 0.4;
padFact     = 5;
scalef      = 1;
maskDiam    = 1550./scalef; %1650            % diameter mask size in um

pixsize1 = 1.44276;             % um            % um           % pixelsize in um
pixsize2 = 1.44978;
pixsize  = pixsize2;

% distance as in the manuscript of Seaberg
R1 = 478e-3;        % in meters
R2 = 1925e-3;       % in meters

magnif = 11;

per_diag = sqrt(2)*magnif*PitchG;


saveName = 'E:\LCLS\processing\Grating_1s\Checkboard_Plate.mat';

if undersampl == 0, scalef = 1;end
%% ===================================================================================
%                                 calculate mask
% ====================================================================================
%return;
[m, n, nImages] = size(stsample);

file1 = mean(stsample,3);
if ~exist('Center','var'),     Center = [];        end
[smmask, bigmask,Center] = maskbuilder(file1,maskDiam/pixsize,Center);
[JH, JV] = meshgrid(1:size(bigmask,2),1:(size(bigmask,1)));
file1 = file1(min(JV(bigmask==1)):max(JV(bigmask==1)) , min(JH(bigmask==1)):max(JH(bigmask==1)));

%% ===================== Fourier transform ===============================%%

ftstack = zeros(m,n,nImages);
for k = 1:nImages
    ftstack(:,:,k) = fft2(stsample(:,:,k));
end;
figure(18)
imagesc(fftshift(sum(log(abs(ftstack)),3)))
colormap gray
%% find centers peak 1 and 2
piece1 = sum(abs(ftstack(1:round(2.5*per_diag/pixsize),1:round(2.5*per_diag/pixsize),:)),3);
piece2 = sum(abs(ftstack(1:round(2.5*per_diag/pixsize),end-round(2.5*per_diag/pixsize):end,:)),3);
piece1(1:10,:) = piece1(end,end);  piece1(:,1:10) = piece1(end,end);
piece2(1:10,:) = piece2(end,1);piece2(:,end-10:end) = piece2(end,1);

[m1,n1] = size(piece1);
[X,Y] = meshgrid(1:n1,1:m1);
X = X - round(m1/2);     Y = Y - round(n1/2);
peaksz = PitchG*magnif/pixsize*sqrt(2)/5;
elpeak = (X.^2 + Y.^2)./(2*peaksz);
elpeak       = exp(-elpeak);

h1 = conv2(piece1,elpeak,'same');       h2 = conv2(piece2,elpeak,'same');

[xpeak1, ypeak1, max_f1] = findpeak(h1,1);
[xpeak2, ypeak2, max_f2] = findpeak(rot90(h2,1),1);

%% ===========================    Find gradient =========================
imk = stsample(:,:,1);

blackman2d = repmat(blackman(m),1,n).* repmat(blackman(n)',m,1);

ftimk = fft2(imk.*blackman2d,m*5,n*5);
ftimk(round(ypeak1*1.5*5):end,:) = 0;    ftimk(:,round(xpeak1*1.5*5):end) = 0;
ftimk(1:round(ypeak1*0.5*5),:) = 0;      ftimk(:,1:round(xpeak1*0.5*5)) = 0;
ftimk = circshift(circshift(ftimk,-round(ypeak1*5),1),-round(xpeak1*5),2);


imkinv = angle(ifft2(ftimk));
imkinv1 = imkinv(1:m,1:n);

ftimk = fft2(rot90(imk,1).*blackman2d',n*5,m*5);
ftimk(round(ypeak2*1.5*5):end,:) = 0;         ftimk(:,round(xpeak2*1.5*5):end) = 0;
ftimk(1:round(ypeak2*0.5*5),:) = 0;           ftimk(:,1:round(xpeak2*0.5*5)) = 0;
ftimk = circshift(circshift(ftimk,-round(ypeak2*5),1),-round(xpeak2*5),2);


imkinv = angle(ifft2(ftimk));
imkinv2 = rot90(imkinv(1:n,1:m),-1);

Phase = cat(3,imkinv1,imkinv2);
Phase = Phase(min(JV(bigmask==1)):max(JV(bigmask==1)) , min(JH(bigmask==1)):max(JH(bigmask==1)),:);

%% ========================================================================

phaseMap = Phase(:,:,1).*smmask;
[m1,n1] = size(phaseMap);

p1 = fliplr(flipud(unwrap(unwrap(fliplr(flipud(phaseMap(1:floor(m1/2),1:floor(n1/2)))),[],2 )))); %#ok<*FLUDLR>
p2 = flipud(unwrap(unwrap(flipud(phaseMap(1:floor(m1/2),floor(n1/2)+1:end)),[],2 )));
p3 = fliplr(unwrap(unwrap(fliplr(phaseMap(floor(m1/2)+1:m1,1:floor(n1/2))),[],2 )));
p4 = unwrap(unwrap(phaseMap(floor(m1/2)+1:m1,floor(n1/2)+1:n1),[],2 ));

dv1 = [p1 p2 ; p3 p4];

%%%   second direction
phaseMap = Phase(:,:,2).*smmask;
[m1,n1] = size(phaseMap);

p1 = fliplr(flipud(unwrap(unwrap(fliplr(flipud(phaseMap(1:floor(m1/2),1:floor(n1/2)))),[],2 )))); %#ok<*FLUDLR>
p2 = flipud(unwrap(unwrap(flipud(phaseMap(1:floor(m1/2),floor(n1/2)+1:end)),[],2 )));
p3 = fliplr(unwrap(unwrap(fliplr(phaseMap(floor(m1/2)+1:m1,1:floor(n1/2))),[],2 )));
p4 = unwrap(unwrap(phaseMap(floor(m1/2)+1:m1,floor(n1/2)+1:n1),[],2 ));

dv2 = [p1 p2 ; p3 p4];

[reconst1, reconst2] = gradient_error(dv1,dv1,smmask);
dv1p = dv1 - reconst1 - reconst2;
dv1p = dv1p - mean(dv1p(:));
[reconst1, reconst2] = gradient_error(dv2,dv2,smmask);
dv2p = dv2 - reconst1 - reconst2;
dv2p = dv2p - mean(dv2p(:));


%% ========================= plots  ==================================%%
figure(1)
subplot(1,2,1);imagesc(dv1p.*smmask,[-1 1].*pi)
colorbar('SouthOutSide')
subplot(1,2,2);imagesc(dv2p.*smmask,[-1 1].*pi)
colorbar('SouthOutSide')
colormap('winter')



figure(8)
subplot(1,2,1)
imagesc(imkinv1)
subplot(1,2,2)
imagesc(imkinv2)

%% =============================  project on the axis ======================

yper1 = m/ypeak1;       xper1 = n/xpeak1; % use period
angleRot = -atan(ypeak1/xpeak1);

dv1pd = dv1p ./(2*pi) .* ((m/(ypeak1-1))^2 + (n/(xpeak1-1))^2 ).^(1/2)/2 .*pixsize*1e-6 ; % scale with the period of the peak in Fourier
dv2pd = dv2p ./(2*pi) .* ((m/(ypeak1-1))^2 + (n/(xpeak1-1))^2 ).^(1/2)/2 .*pixsize*1e-6 ; % also why the 2 here

gradH = (-cos(angleRot)*dv1pd - sin(angleRot) *dv2pd);  % rotation
gradV = (sin(angleRot) *dv1pd - cos(angleRot) *dv2pd);
%% =============================  scale the phase ======================
% This is the formula from seaberg
gradH = gradH.*R1./R2./(R2 - R1)./2;  % wavefront gradient and not phase gradient
gradV = gradV.*R1./R2./(R2 - R1)./2;  % why the two?

if ~isempty(saveName),  save(saveName,'gradH','gradV','file1','smmask','Center');end;
return;

figure(23)
subplot(2,3,1)
imagesc( dv1.*smmask,[-1 1].*dispscale )
set(gca,'DataAspectRatio',[1 1 1]')
set(gca,'YTick',0);set(gca,'XTick',0);  set(gca,'TickLength',[0 0]);
t = colorbar('SouthOutside','peer',gca);
set(get(t,'ylabel'),'String', 'a.u','FontName','Times','FontSize',14);
title('GradV (horiz direc)')
% scalebar('Location','southeast','Colour',barcolor,'Bold','True','Unit','mm','ScaleLength',n2)
title('GradV (vert direc)')

subplot(2,3,2)
imagesc(dv2.*smmask,[-1 1].*dispscale )
set(gca,'DataAspectRatio',[1 1 1]')
set(gca,'YTick',0);set(gca,'XTick',0);  set(gca,'TickLength',[0 0]);
t = colorbar('SouthOutside','peer',gca);
set(get(t,'ylabel'),'String', 'a.u','FontName','Times','FontSize',14);
title('GradV (horiz direc)')
% scalebar('Location','southeast','Colour',barcolor,'Bold','True','Unit','mm','ScaleLength',n2)



subplot(2,3,4)
imagesc(dh1.*smmask,[-1 1.*dispscale ])% + hdir*pixreach.*pixsize/step)
set(gca,'DataAspectRatio',[1 1 1]')
set(gca,'YTick',0);set(gca,'XTick',0);  set(gca,'TickLength',[0 0]);
t = colorbar('SouthOutside','peer',gca);
set(get(t,'ylabel'),'String', 'a.u','FontName','Times','FontSize',14);
title('GradH (vert direc)')
% scalebar('Location','southeast','Colour',barcolor,'Bold','True','Unit','mm','ScaleLength',n2)


subplot(2,3,5)
imagesc(dh2.*smmask,[-1 1].*dispscale )
set(gca,'DataAspectRatio',[1 1 1]')
set(gca,'YTick',0);set(gca,'XTick',0);  set(gca,'TickLength',[0 0]);
t = colorbar('SouthOutside','peer',gca);
set(get(t,'ylabel'),'String', 'a.u','FontName','Times','FontSize',14);
title('GradH (horiz direc)')
% scalebar('Location','southeast','Colour',barcolor,'Bold','True','Unit','mm','ScaleLength',n2)


subplot(2,3,3)
%plot(mean(GradV(:,:,2)))% + vdir*pixreach.*pixsize/step,2))
z1 = frankotchellappa(dv1.*smmask,dv2.*smmask);
z1 = z1 - mean(z1(smmask(:)));
imagesc(z1.*smmask)

subplot(2,3,6)
%plot(mean(GradH(:,:,1)))%+ hdir*pixreach.*pixsize/step,1))
z2 = frankotchellappa(dh1.*smmask,dh2.*smmask);
z2 = z2 - mean(z2(smmask(:)));
imagesc(z2.*smmask)
colormap(othercolor('RdYlGn10'))
colormap gray

z3 = frankotchellappa(z2.*smmask,z1.*smmask);

% set(23,'PaperPositionMode','auto','PaperType','A3','PaperOrientation','landscape')
% print(23,'-dsvg',fullfile(['/users/berujon/Report+Publi/SpeckleMetrol/figures/curMirr.svg']));

