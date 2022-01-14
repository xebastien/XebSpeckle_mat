% ===============================================================
pixreach = 6; % number of interpixel distance to correlate
scanDimOrd = 1; % )How was made the mesh scan: 1st dim is vert or horizontal ? boolean 1 or 0
ROIwidth =  18;


subdiv = 19;
select_ROI = 1;         % 0 = full field, 1 = ROI selection, 2 = manual below 3-auto    4-reuse from 1
%[top down left right]  = [mindim1 maxdim1 mindim2 maxdim2] when imagesc
ROI1 = [470 980 235 727].*2;% LCLS first detector
%ROI1 = [675 775 440 522].*2./scalef;% LCLS first detector
%ROI1 = [875 900 463 500].*2./scalef;% LCLS first detector
%ROI2 = [240 830 220 800].*2./scalef;% LCLS second detector
MaxDisp = 11;        % define the dynamic of the method - high = high dynamic by sensitivity may reduce

clear delaymat
dispscale = 1/pixreach;
%% ========================================================================
%                       process vertical direction                        %
% =========================================================================
offPR = floor(pixreach/2);

tic
stsampleV1 = stsampleV(ROI1(1) +round(ROIwidth/2)+1:ROI1(2)-round(ROIwidth/2),ROI1(3):ROI1(4),:);
% offPR:ROI1(2),ROI1(3) + offPR:ROI1(4) - offPR,:);
%stsample1 = stref(ROI2(1):ROI2(2),ROI2(3):ROI2(4),:);
[m1,n1,r1] = size(stsampleV1);
%for k = 1:size(stsample1,1), stsample1(k,:,:) = medfilt2(squeeze(stsample1(k,:,:)),[3 3]);end;
% for k = 1:size(stsample1,1), stsample(:,:,k) = (single(stsample(:,:,k)) - mean(mean(single(stsample(:,:,k)))))./std(single(reshape(stsample(:,:,k),1,[])));end;
for k = 1:size(stsampleV1,1), stsampleV1(k,:,:) = myerosion(squeeze(stsampleV1(k,:,:)),1,5,'silent');end;

ROIcent = 1 + round(ROIwidth/2) + 1 :1: n1 - round(ROIwidth/2)-1;
delaymatV = zeros(m1-pixreach,2,length(ROIcent));
%return;
for kROI =1 :1: length(ROIcent)
    ROI = round([1 m1 (ROIcent(kROI)- round(ROIwidth/2)) (ROIcent(kROI)+ round(ROIwidth/2))]);
    % make a stack of sample images
    stack1 = stsampleV1(ROI(1):ROI(2),ROI(3):ROI(4),:);
    
    clear delay coeff pic      
    delay = zeros(m1-pixreach,2);     coeff = zeros(m1-pixreach,1);     pic = zeros(m1-pixreach,5);

    parfor pp = 1 :1: m1-pixreach
       [delay(pp,:),coeff(pp),pic(pp,:)] = SliceDelayXdefl(squeeze(double(stack1(pp,:,:))), squeeze(double(stack1(pp+pixreach,:,:))),MaxDisp);
       %progmeter(pp/(m1-nstep1));
    end;
    delaymatV(:,:,kROI) = delay;
    disp(['Loop  ' num2str(kROI) ' done']);
end;


dmVV = squeeze(delaymatV(:,2,:));             dmVH = squeeze(delaymatV(:,1,:));
dmVV = dmVV(:,offPR:end-offPR);               dmVH = dmVH(:,offPR:end-offPR);
%% ========================================================================
%                       process horizontal direction                       %
% =========================================================================


stsampleH1 = stsampleH(ROI1(1):ROI1(2),ROI1(3) + round(ROIwidth/2)+1:ROI1(4)-round(ROIwidth/2),:);
%ROI1(1) + offPR:ROI1(2) - offPR,ROI1(3) + offPR:ROI1(4),:);
stsampleH1 = permute(stsampleH1,[2 1 3]);
[m1,n1,r1] = size(stsampleH1);

for k = 1:size(stsampleH1,1), stsampleH1(k,:,:) = myerosion(squeeze(stsampleH1(k,:,:)),1,5,'silent');end;

ROIcent = 1 + round(ROIwidth/2) + 1 :1: n1 - round(ROIwidth/2)-1;
delaymatH = zeros(m1-pixreach,2,length(ROIcent));

for kROI =1 :1: length(ROIcent)
    ROI = round([1 m1 (ROIcent(kROI)- round(ROIwidth/2)) (ROIcent(kROI)+ round(ROIwidth/2))]);
    % make a stack of sample images
    stack1 = stsampleH1(ROI(1):ROI(2),ROI(3):ROI(4),:);
    
    clear delay coeff pic      
    delay = zeros(m1-pixreach,2);     coeff = zeros(m1-pixreach,1);     pic = zeros(m1-pixreach,5);

    parfor pp = 1 :1: m1-pixreach
       [delay(pp,:),coeff(pp),pic(pp,:)] = SliceDelayXdefl(squeeze(double(stack1(pp,:,:))), squeeze(double(stack1(pp+pixreach,:,:))),MaxDisp);
       %progmeter(pp/(m1-nstep1));
    end;
    delaymatH(:,:,kROI) = delay;
    disp(['Loop  ' num2str(kROI) ' done']);
end;
   
dmHH = squeeze(delaymatH(:,2,:))';            dmHV = squeeze(delaymatH(:,1,:))';
dmHH = dmHH(offPR:end-offPR,:);               dmHV = dmHV(offPR:end-offPR,:);
% %% #########################################################################

msmin = min([size(dmVV);size(dmHH)]);

dmHH = dmHH(1:msmin(1),1:msmin(2));     dmVV = dmVV(1:msmin(1),1:msmin(2));
dmHV = dmHV(1:msmin(1),1:msmin(2));     dmHV = dmHV(1:msmin(1),1:msmin(2));


dmHH(isnan(dmHH)) = mean(dmHH(~isnan(dmHH)));
dmVH(isnan(dmVH)) = mean(dmVH(~isnan(dmVH)));
dmHV(isnan(dmHV)) = mean(dmHV(~isnan(dmHV)));
dmVV(isnan(dmVV)) = mean(dmVV(~isnan(dmVV)));

dmHH = myerosion(dmHH,3,3);     dmVH = myerosion(dmVH,3,3);
dmHV = myerosion(dmHV,3,3);     dmVV = myerosion(dmVV,3,3);

dmHH = dmHH./pixreach;         dmHV = dmHV./pixreach;
dmVV = dmVV./pixreach;         dmVH = dmVH./pixreach;



figure(2)
subplot(2,2,1)
imagesc(dmVV,([-2 2]./3./pixreach + median(dmVV(:))))
title('dmVV')
colorbar

subplot(2,2,2)
imagesc(dmVH,[-2 2]./8./pixreach + median(dmVH(:)))
title('dmVH')
colorbar

subplot(2,2,3)
imagesc(dmHH,[-2 2]./3./pixreach + median(dmHH(:)))
title('dmHH')
colorbar

subplot(2,2,4)
imagesc(dmHV,[-2 2]./8./pixreach + median(dmHV(:)))
title('dmHV')
colorbar

colormap gray


%% ======================================================================%%
barcolor = [1 1 1].*0;
pixsize = 1.44276;             % um            % um           % pixelsize in um
maskDiam = 1270;            % diameter mask size in um

dispscale = 0.5;


file1 = mean(stsampleV1(1+ offPR:end- offPR,round(ROIwidth/2)+1+ offPR:end-round(ROIwidth/2)- offPR,:),3);
file1 = file1(1:msmin(1),1:msmin(2));  
[smmask, bigmask,Center] = maskbuilder(file1,maskDiam/pixsize);

[JH, JV] = meshgrid(1:size(bigmask,2),1:(size(bigmask,1)));
dv1 = dmVV(min(JV(bigmask==1)):max(JV(bigmask==1)) , min(JH(bigmask==1)):max(JH(bigmask==1)));
dv2 = dmVH(min(JV(bigmask==1)):max(JV(bigmask==1)) , min(JH(bigmask==1)):max(JH(bigmask==1)));

dh1 = dmHV(min(JV(bigmask==1)):max(JV(bigmask==1)) , min(JH(bigmask==1)):max(JH(bigmask==1)));
dh2 = dmHH(min(JV(bigmask==1)):max(JV(bigmask==1)) , min(JH(bigmask==1)):max(JH(bigmask==1)));

if mean(dv1(smmask)) < 0 ,dv1 = -dv1;dv2 = -dv2;end;
if mean(dh2(smmask)) < 0 ,dh1 = -dh1;dh2 = -dh2;end;
% dv1 = myerosion(dv1,3,3);       dv2 = myerosion(dv2,3,3); dh1 =
% myerosion(dh1,3,3);       dh2 = myerosion(dh2,3,3);

% dv1 = dv1 - mean(dv1(smmask(:)));           dv2 = dv2 - mean(dv2(smmask(:)));
% dh1 = dh1 - mean(dh1(smmask(:)));           dh2 = dh2 - mean(dh2(smmask(:)));

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

set(23,'PaperPositionMode','auto','PaperType','A3','PaperOrientation','landscape')
print(23,'-dsvg',fullfile(['/users/berujon/Report+Publi/SpeckleMetrol/figures/curMirr.svg']));



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
