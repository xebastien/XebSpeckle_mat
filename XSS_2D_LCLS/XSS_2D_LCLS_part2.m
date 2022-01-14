% ===============================================================
pixreach.V = 5; % number of interpixel distance to correlate
pixreach.H = 5;
scanDimOrd = 1; % )How was made the mesh scan: 1st dim is vert or horizontal ? boolean 1 or 0


if undersampl == 0, scalef = 1;end
subdiv = 19;
select_ROI = 1;         % 0 = full field, 1 = ROI selection, 2 = manual below 3-auto    4-reuse from 1
%[top down left right]  = [mindim1 maxdim1 mindim2 maxdim2] when imagesc
ROI1 = [475 975 240 722].*2./scalef;% LCLS first detector

%ROI1 = [240 830 220 800].*2./scalef;% LCLS second detector
%ROI1 = [520 1634 490 1605];
%ROI1 = [800 1750 400 1700];% LCLS second detector grating

ROI1 = [950 1950 480 1444];% LCLS first detector speckle / no plate
ROI1 = [320 1934 290 1805];% LCLS second detector speckle / with and without plate

dispscale = 1/pixreach.V;

saveEnvironment = '/data/bm05/imaging/seb/LCLS/processing/XSS2D/'; 
enviname = 'plateSpeckleDet1';
%% ===================================================================================
%                                   Build cells
% ====================================================================================
%return;
[m1,n1] = size(files2(1).data);
zonebounds = round(linspace(1,m1-pixreach.V,subdiv+1));

GradHc = cell(1,subdiv);    GradVc = cell(1,subdiv);
GradD1c = cell(1,subdiv);   GradD2c = cell(1,subdiv);
piece_stack = cell(1,subdiv);
%% ======================================================================%%
tic
stsample1 = stsample(ROI1(1):ROI1(2),ROI1(3):ROI1(4),:);
stsample1 = stref(ROI1(1):ROI1(2),ROI1(3):ROI1(4),:);
%
%stsample1 = stref(ROI2(1):ROI2(2),ROI2(3):ROI2(4),:);
[m1,n1,r1] = size(stsample1);
%clear edgeSize;
%return;
if sqrt(nImages) ~= round(sqrt(nImages)) && sqrt(nImages-1) == round(sqrt(nImages-1)),
    nImages = nImages - 1;
end;
zonebounds = round(linspace(1,m1-pixreach.V,subdiv+1));
for k = 1 : subdiv
        if k==1, zone = zonebounds(k):1:zonebounds(k+1);
        else     zone = zonebounds(k) + 1:1:zonebounds(k+1);
        end;
        zone = [zone (zonebounds(k+1)+1 : 1 :zonebounds(k+1) + pixreach.V) ];
        piece_stack{k} = zeros(length(zone) , n1, nImages );
        for pp = 1 : nImages, 
            piece_stack{k}(:,:,pp) = stsample1(zone,:,pp); 
        end;
end;%#ok<PFBNS>


IVc = zeros(subdiv,2);      IHc = zeros(subdiv,2);
ID1 = zeros(subdiv,2);      ID2 = zeros(subdiv,2);
parfor k = 1 : subdiv    
   IVc(k,:) = pixdelaycorr_findI0(piece_stack{k}(1:end-pixreach.V,1:end-pixreach.H,:),piece_stack{k}(pixreach.V+1:end,1:end-pixreach.H,:));
   IHc(k,:) = pixdelaycorr_findI0(piece_stack{k}(1:end-pixreach.V,1:end-pixreach.H,:),piece_stack{k}(1:end-pixreach.V,pixreach.H+1:end,:));
%    ID1(k,:) = pixdelaycorr_findI0(piece_stack{k}(1:end-pixreach.V,1:end-pixreach.H,:),piece_stack{k}(pixreach.V+1:end,pixreach.H+1:end,:));     
%    ID2(k,:) = pixdelaycorr_findI0(piece_stack{k}(1:end-pixreach.V,pixreach.H+1:end,:),piece_stack{k}(pixreach.V+1:end,1:end-pixreach.H,:));
   % progmeter(k/(subdiv))
end;
IVc = round(mean(IVc));IHc = round(mean(IHc));
ID1 = round(mean(ID1));ID2  = round(mean(ID2));
disp(['IVc  ' num2str(IVc) '   IHc  ' num2str(IHc)]);
disp('-----------------')
    % -------- correlations ---------%
%if parpool('size')>0
parfor k = 1 : subdiv    
   GradVc{k} = pixdelaycorr2(piece_stack{k}(1:end-pixreach.V,1:end-pixreach.H,:),piece_stack{k}(pixreach.V+1:end,1:end-pixreach.H,:),IVc);
   GradHc{k} = pixdelaycorr2(piece_stack{k}(1:end-pixreach.V,1:end-pixreach.H,:),piece_stack{k}(1:end-pixreach.V,pixreach.H+1:end,:),IHc);
   GradD1c{k} = pixdelaycorr2(piece_stack{k}(1:end-pixreach.V,1:end-pixreach.H,:),piece_stack{k}(pixreach.V+1:end,pixreach.H+1:end,:),ID1);     
%    GradD2c{k} = pixdelaycorr2(piece_stack{k}(1:end-pixreach.V,pixreach.H+1:end,:),piece_stack{k}(pixreach.V+1:end,1:end-pixreach.H,:),ID2);
   progmeter(k/(subdiv))
end;
toc

GradV = cat(1,GradVc{:});
GradH = cat(1,GradHc{:});
GradD1 = cat(1,GradD1c{:}); 
GradD2 = cat(1,GradD2c{:}); 

% --------------- integration to get the wavefront gradient --------------%
GradVs = GradV;         GradHs = GradH;
GradD1s = GradD1;       GradD2s = GradD2;
GradH1 = GradH;         GradV1 = GradV;% backups   
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
maskDiam = 1550./scalef;            % diameter mask size in um

dispscale = 0.5;

file1 = mean(stsample1(1:end-pixreach.V,1:end-pixreach.H,:),3);
[smmask, bigmask,Center] = maskbuilder(file1,maskDiam/pixsize2);

[JH, JV] = meshgrid(1:size(bigmask,2),1:(size(bigmask,1)));
dv1 = GradV(min(JV(bigmask==1)):max(JV(bigmask==1)) , min(JH(bigmask==1)):max(JH(bigmask==1)),1);
dv2 = GradV(min(JV(bigmask==1)):max(JV(bigmask==1)) , min(JH(bigmask==1)):max(JH(bigmask==1)),2);

dh1 = GradH(min(JV(bigmask==1)):max(JV(bigmask==1)) , min(JH(bigmask==1)):max(JH(bigmask==1)),1);
dh2 = GradH(min(JV(bigmask==1)):max(JV(bigmask==1)) , min(JH(bigmask==1)):max(JH(bigmask==1)),2);

% dv1 = myerosion(dv1,3,3);       dv2 = myerosion(dv2,3,3);
% dh1 = myerosion(dh1,3,3);       dh2 = myerosion(dh2,3,3);

dv1 = dv1 - mean(dv1(smmask(:)));           dv2 = dv2 - mean(dv2(smmask(:)));
dh1 = dh1 - mean(dh1(smmask(:)));           dh2 = dh2 - mean(dh2(smmask(:)));



if ~isempty(saveEnvironment)
    save(fullfile(saveEnvironment,[enviname 'mat']),'dv1','dv2','file1','dh1','dh2','maskDiam','smmask','Center');
end


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
z1 = frankotchellappa(dv2.*smmask,-dv1.*smmask);
z1 = z1 - mean(z1(smmask(:)));
imagesc(z1.*smmask)

subplot(2,3,6)
%plot(mean(GradH(:,:,1)))%+ hdir*pixreach.*pixsize/step,1))
z2 = frankotchellappa(dh2.*smmask,-dh1.*smmask);
z2 = z2 - mean(z2(smmask(:)));
imagesc(z2.*smmask)
colormap(othercolor('RdYlGn10'))
colormap gray

z3 = frankotchellappa(z2.*smmask,z1.*smmask);

% set(23,'PaperPositionMode','auto','PaperType','A3','PaperOrientation','landscape')
% print(23,'-dsvg',fullfile(['/users/berujon/Report+Publi/SpeckleMetrol/figures/curMirr.svg']));



return;
figure(24)
subplot(2,3,1)
imagesc(GradH(:,:,2))
title(['Avg  H = ' num2str(mean2(GradH1(:,:,1)))]);
colorbar
subplot(2,3,2)
imagesc(GradD1(:,:,1))
title(['Avg H-D1= ' num2str(mean2(GradD1(:,:,1)))]);
colorbar
subplot(2,3,3)
imagesc(-GradD2(:,:,1))
title(['Avg H-D2= ' num2str(mean2(GradD2(:,:,1)))]);
colorbar

subplot(2,3,4)
imagesc(GradV(:,:,1))
title(['Avg V = ' num2str(mean2(GradV1(:,:,2)))]);
colorbar
subplot(2,3,5)
imagesc(GradD1(:,:,2))
title(['Avg V-D1 = ' num2str(mean2(GradD1(:,:,2)))]);
colorbar
subplot(2,3,6)
imagesc(GradD2(:,:,2))
title(['Avg V-D2= ' num2str(mean2(GradD2(:,:,2)))]);
colorbar



%GradV2 = intgrad2(GradV2(:,:,1),GradV2(:,:,2),1,1);
%GradH2 = intgrad2(GradH2(:,:,1),GradH2(:,:,2),1,1);
return;
