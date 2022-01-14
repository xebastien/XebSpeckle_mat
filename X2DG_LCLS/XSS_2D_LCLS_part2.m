% ===============================================================
pixreach.V = 3; % number of interpixel distance to correlate
pixreach.H = 3;
scanDimOrd = 1; % )How was made the mesh scan: 1st dim is vert or horizontal ? boolean 1 or 0
PitchG = 4; StepS = 0.4;
padFact =5;
scalef = 1;
maskDiam = 1650./scalef;             % diameter mask size in um

pixsize1 = 1.44276;             % um            % um           % pixelsize in um
pixsize2 = 1.44978;
pixsize = pixsize2;


if undersampl == 0, scalef = 1;end
subdiv = 59;        % for parallel processing

%[top down left right]  = [mindim1 maxdim1 mindim2 maxdim2] when imagesc
ROI1 = [480 1641 320 1470];         % LCLS second detector
%ROI1 = [450 1691 280 1520];        % LCLS second detector scans69-95
ROI1 = [480 1641 320+93 1470+160];  % LCLS second detector 
 ROI1 = [480 1641 320 1470];  % LCLS second detector               Phase Plate or not +  Diamond Checkboard 4um
% for debug (inside of the beam/ no dark in aimge)
%ROI1 = [780 1241 720+93 1170+93];   % LCLS second detector


dispscale = 1/pixreach.V;
%% ===================================================================================
%                                   Build cells
% ====================================================================================
%return;
[m1,n1] = size(files{1}(1).data);
zonebounds = round(linspace(1,m1-pixreach.V,subdiv+1));

GradHc = cell(1,subdiv);    GradVc = cell(1,subdiv);
GradD1c = cell(1,subdiv);   GradD2c = cell(1,subdiv);
piece_stack = cell(1,subdiv);
%% ======================================================================%%
tic
stsample1 = stsample(ROI1(1):ROI1(2),ROI1(3):ROI1(4),:);
%
%stsample1 = stref(ROI2(1):ROI2(2),ROI2(3):ROI2(4),:);
[m1,n1,r1] = size(stsample1);
%clear edgeSize;

if sqrt(nImages) ~= round(sqrt(nImages)) && sqrt(nImages-1) == round(sqrt(nImages-1)),
    nImages = nImages - 1;
end;
zonebounds = round(linspace(1,m1,subdiv+1));
for k = 1 : subdiv
        if k==1, zone = zonebounds(k):1:zonebounds(k+1);
        else     zone = zonebounds(k) + 1:1:zonebounds(k+1);
        end;
        zone = [zone (zonebounds(k+1)+1 : 1 :zonebounds(k+1) ) ];
        piece_stack{k} = zeros(length(zone) , n1, nImages );
        for pp = 1 : nImages, 
            piece_stack{k}(:,:,pp) = stsample1(zone,:,pp); 
        end;
end;%#ok<PFBNS>

nStep = sqrt(size(stsample,3));
PeriodScan = StepS*nStep/PitchG;

%if parpool('size')>0
Amp = cell(1,subdiv); Phase = cell(1,subdiv); Scat = cell(1,subdiv); 
for k = 1 : subdiv    
   [Amp{k},Phase{k},Scat{k}] = psiphase2D(piece_stack{k},PeriodScan,padFact);
end
toc

Amp  = cat(1,Amp{:});Phase  = cat(1,Phase{:});Scat  = cat(1,Scat{:});
%GradH  = cat(1,GradHc{:});
% GradD1 = cat(1,GradD1c{:}); 
% GradD2 = cat(1,GradD2c{:}); 

% --------------- integration to get the wavefront gradient --------------%

AmpS = Amp ;PhaseS = Phase; ScatS = Scat;

% calculate mask

file1 = mean(stsample1,3);

if ~exist('Center','var'),     Center = [];        end
[smmask, bigmask,Center] = maskbuilder(file1,maskDiam/pixsize,Center);


[JH, JV] = meshgrid(1:size(bigmask,2),1:(size(bigmask,1)));
file1 = file1(min(JV(bigmask==1)):max(JV(bigmask==1)) , min(JH(bigmask==1)):max(JH(bigmask==1)));
Phase = Phase(min(JV(bigmask==1)):max(JV(bigmask==1)) , min(JH(bigmask==1)):max(JH(bigmask==1)),:);



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

save('/data/bm05/imaging/seb/LCLS/processing/Grating/Checkboard_Plate.mat','dv1','AmpS','ScatS','dv2','file1','smmask','Center','bigmask')
% 


[reconst1, reconst2] = gradient_error(dv1,dv1,smmask);
dv1p = dv1 - reconst1 - reconst2;
dv1p = dv1p - mean(dv1p(:));
[reconst1, reconst2] = gradient_error(dv2,dv2,smmask);
dv2p = dv2 - reconst1 - reconst2;
dv2p = dv2p - mean(dv2p(:));




figure(1)
subplot(1,2,1);imagesc(dv1p.*smmask,[-1 1].*1)
subplot(1,2,2);imagesc(dv2p.*smmask,[-1 1].*1)
return;
clear GradV GradH

if scanDimOrd,
    GradV(:,:,1) = GradV1(:,:,2) ./ pixreach.V;   GradH(:,:,1) = GradH1(:,:,2) ./ pixreach.H;
    GradV(:,:,2) = GradV1(:,:,1) ./ pixreach.V;   GradH(:,:,2) = GradH1(:,:,1) ./ pixreach.H;
else
    GradV(:,:,1) = GradH1(:,:,1) ./ pixreach.H;   GradH(:,:,1) = GradV1(:,:,2) ./ pixreach.V;
    GradV(:,:,2) = GradH1(:,:,2) ./ pixreach.H;   GradH(:,:,2) = GradV1(:,:,2) ./ pixreach.V;
end;
    
clear GradVc GradHc GradD1c GradD2c

%% ======================================================================%%
n2 = 172;
barcolor = [1 1 1].*0;
pixsize = 1.44276;             % um            % um           % pixelsize in um
pixsize2 = 1.44978;

dispscale = 0.5;

file1 = mean(stsample1(1:end-pixreach.V,1:end-pixreach.H,:),3);
[smmask, bigmask,Center] = maskbuilder(file1,maskDiam/pixsize);

[JH, JV] = meshgrid(1:size(bigmask,2),1:(size(bigmask,1)));
dv1 = GradV(min(JV(bigmask==1)):max(JV(bigmask==1)) , min(JH(bigmask==1)):max(JH(bigmask==1)),1);
dv2 = GradV(min(JV(bigmask==1)):max(JV(bigmask==1)) , min(JH(bigmask==1)):max(JH(bigmask==1)),2);

dh1 = GradH(min(JV(bigmask==1)):max(JV(bigmask==1)) , min(JH(bigmask==1)):max(JH(bigmask==1)),1);
dh2 = GradH(min(JV(bigmask==1)):max(JV(bigmask==1)) , min(JH(bigmask==1)):max(JH(bigmask==1)),2);

% dv1 = myerosion(dv1,3,3);       dv2 = myerosion(dv2,3,3);
% dh1 = myerosion(dh1,3,3);       dh2 = myerosion(dh2,3,3);

dv1 = dv1 - mean(dv1(smmask(:)));           dv2 = dv2 - mean(dv2(smmask(:)));
dh1 = dh1 - mean(dh1(smmask(:)));           dh2 = dh2 - mean(dh2(smmask(:)));

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

return;
GradD1(:,:,2) = GradD1(:,:,2) - mean(reshape(GradD1(:,:,2),1,[]));
GradD1(:,:,2) = GradD1(:,:,2) - mean(reshape(GradD1(:,:,2),1,[]));

GradD2(:,:,2) = GradD2(:,:,2) - mean(reshape(GradD2(:,:,2),1,[]));
GradD2(:,:,2) = GradD2(:,:,2) - mean(reshape(GradD2(:,:,2),1,[]));

figure(24)
subplot(2,3,1)
imagesc(GradH(:,:,2).*bigmask,[-1 1].*dispscale.*10 )
title(['Avg  H = ' num2str(mean2(GradH1(:,:,1)))]);
colorbar
subplot(2,3,2)
imagesc(GradD1(:,:,1).*bigmask,[-1 1].*dispscale.*10 )
title(['Avg H-D1= ' num2str(mean2(GradD1(:,:,1)))]);
colorbar
subplot(2,3,3)
imagesc(GradD2(:,:,1).*bigmask,[-1 1].*dispscale.*10 )
title(['Avg H-D2= ' num2str(mean2(GradD2(:,:,1)))]);
colorbar

subplot(2,3,4)
imagesc(GradV(:,:,1).*bigmask,[-1 1].*dispscale )
title(['Avg V = ' num2str(mean2(GradV1(:,:,2)))]);
colorbar
subplot(2,3,5)
imagesc(GradD1(:,:,2).*bigmask,[-1 1].*dispscale)
title(['Avg V-D1 = ' num2str(mean2(GradD1(:,:,2)))]);
colorbar
subplot(2,3,6)
imagesc(GradD2(:,:,2).*bigmask,[-1 1].*dispscale)
title(['Avg V-D2= ' num2str(mean2(GradD2(:,:,2)))]);
colorbar

return;

%GradV2 = intgrad2(GradV2(:,:,1),GradV2(:,:,2),1,1);
%GradH2 = intgrad2(GradH2(:,:,1),GradH2(:,:,2),1,1);
return;
