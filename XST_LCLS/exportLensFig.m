 CM1 = cbrewer('div','RdBu',100);% colormap for global phsae graident
CM2 = colormap('gray');colormap(othercolor('PuBu9'));%cbrewer('seq','PuBu',100);%;colormap(othercolor('Reds9'));%cbrewer('div','RdBu',100);%cbrewer('div','RdGy',100);%cbrewer('div','PuOr',100);%
%colormap(othercolor('RdYlBu11'));% colormap for global phase gradient error
CM4 = cbrewer('div','RdGy',100);%colormap(othercolor('PRGn5'));% colormap for phase error and thickness error
CM3 = colormap(othercolor('RdYlBu11'));%colormap(othercolor('PuBu9')); %colormap('gray');
barcolor = [1 0 0];

%% ----------------figure wavefront gradient as calculated---------------%
[m,n] = size(dphaseX);
m1 = (1:1:m).*pixsize./1000;n1 = (1:1:n).*pixsize./1000;

Center2 = Center.*pixsize./1000;
Center2(2) = Center2(2)- mean(n1);
Center2(1) = Center2(1)- mean(m1);
Center2 =  round(Center2.*1000)./1000;
m1 = m1 - mean(m1);n1 = n1 - mean(n1);

figure(1)
colormap(flipud(CM1))
set(1,'units','points','position',[1800,500,800,400])

subplot(1,2,1)
a = mean(dphaseX(:).*1e6);b = std(dphaseX(:).*1e6);
imagesc(n1,m1,dphaseX.*1e6,[a-3*b a+3*b])
set(gca,'DataAspectRatio',[1 1 1]')
%set(gca,'YTick',0);set(gca,'XTick',0);set(gca,'TickLength',[0 0]);
t = colorbar('SouthOutside','peer',gca);
set(get(t,'ylabel'),'String', 'urad','FontName','Times','FontSize',14);
%scalebar('Location','southeast','Colour',barcolor,'Unit','um')
xlabel(['(mm)  - Lens center at ' num2str(fliplr(Center2(2)))]);
ylabel(['(mm)  - Lens center at ' num2str(fliplr(Center2(1)))]);
title('Vertical wavefront gradient')

subplot(1,2,2)
a = mean(dphaseY(:).*1e6);b = std(dphaseY(:).*1e6);
imagesc(n1,m1,dphaseY.*1e6,[a-3*b a+3*b])
set(gca,'DataAspectRatio',[1 1 1]')
%set(gca,'YTick',0);set(gca,'XTick',0);set(gca,'TickLength',[0 0]);
t = colorbar('SouthOutside','peer',gca);
set(get(t,'ylabel'),'String', 'urad','FontName','Times','FontSize',14);
%scalebar('Location','southeast','Colour',barcolor,'Unit','um')
xlabel(['(mm)  - Lens center at ' num2str(fliplr(Center2(2)))]);
ylabel(['(mm)  - Lens center at ' num2str(fliplr(Center2(1)))]);
title('Horizontal wavefront gradient')

set(1,'PaperPositionMode','auto','PaperType','A3','PaperOrientation','landscape')  


%% ----------------figure wavefront gradient with mask -------------------%
[m,n] = size(xg);
m1 = (1:1:m).*pixsize./1000;n1 = (1:1:n).*pixsize./1000;
m1 = m1 - mean(m1);n1 = n1 - mean(n1);

figure(2)
colormap(flipud(CM2))
set(2,'units','points','position',[2000,400,800,400])

subplot(1,2,1)
a = mean(xg(mask1(:)));b = std(xg(mask1(:)));
imagesc(n1,m1,xg.*mask1,[a-3*b a+3*b])
set(gca,'DataAspectRatio',[1 1 1]')
%set(gca,'YTick',0);set(gca,'XTick',0);set(gca,'TickLength',[0 0]);
t = colorbar('SouthOutside','peer',gca);
set(get(t,'ylabel'),'String', 'urad','FontName','Times','FontSize',14);
%scalebar('Location','southeast','Colour',barcolor,'Unit','um')
xlabel(['(mm) - Focal length = ' num2str(raddiff1) ' m']);
ylabel('(mm)');
title('Vertical wavefront gradient error')

subplot(1,2,2)
a = mean(yg(mask1(:)));b = std(yg(mask1(:)));
imagesc(n1,m1,yg.*mask1,[a-3*b a+3*b])
set(gca,'DataAspectRatio',[1 1 1]')
%set(gca,'YTick',0);set(gca,'XTick',0);set(gca,'TickLength',[0 0]);
t = colorbar('SouthOutside','peer',gca);
set(get(t,'ylabel'),'String', 'urad','FontName','Times','FontSize',14);
%scalebar('Location','southeast','Colour',barcolor,'Unit','um')
xlabel(['(mm) - Focal length = ' num2str(raddiff2) ' m']);
ylabel('(mm)');
title('Horizontal wavefront gradient error')

set(2,'PaperPositionMode','auto','PaperType','A3','PaperOrientation','landscape')  
%% ---------------------------figure masks -------------------------------%
figure(3)
set(3,'units','points','position',[1800,200,800,400])

subplot(1,2,1)
imagesc((~mask2).*thephase(:,:,1))
set(gca,'DataAspectRatio',[1 1 1]')
set(gca,'YTick',0);set(gca,'XTick',0);set(gca,'TickLength',[0 0]);
title(['Mask diameter = ' num2str(printmaskSizeV)]);
xlabel(['Center location in pixel ' num2str(fliplr(Center))]);

subplot(1,2,2)
imagesc(d1+mask2.*0.5)
set(gca,'DataAspectRatio',[1 1 1]')
set(gca,'YTick',0);set(gca,'XTick',0);set(gca,'TickLength',[0 0]);
title(['Mask diameter = ' num2str(printmaskSizeV)]);
xlabel(['Center location in pixel ' num2str(fliplr(Center))]);

colormap(CM3)

set(3,'PaperPositionMode','auto','PaperType','A3','PaperOrientation','landscape')  


%% --------------figure wavefront error and thickness error --------------%
%% ----------------figure wavefront gradient with mask -------------------%
[m,n] = size(phi);
m1 = (1:1:m).*pixsize./1000;n1 = (1:1:n).*pixsize./1000;
m1 = m1 - mean(m1);n1 = n1 - mean(n1);

figure(4)
set(4,'units','points','position',[1700,1,1200,400])
set(4,'PaperPositionMode','auto','PaperType','A3','PaperOrientation','landscape') 
colormap(flipud(CM4))

subplot(2,3,1)
a = mean(phi(mask1(:)));b = std(phi(mask1(:)));
imagesc(n1,m1,phi.*mask1.*1e9,[a-3*b a+3*b].*1e9)
set(gca,'DataAspectRatio',[1 1 1]')
%set(gca,'YTick',0);set(gca,'XTick',0);set(gca,'TickLength',[0 0]);
t = colorbar('SouthOutside','peer',gca);
set(get(t,'ylabel'),'String', 'nm','FontName','Times','FontSize',14);
%scalebar('Location','southeast','Colour',barcolor,'Unit','um')
xlabel('(mm)'); ylabel('(mm)');
title(['Wavefront error - RMS = ' num2str(std( phi(mask1(:))*1e9 )) ' nm'] )


subplot(2,3,2)
a = mean(phi(mask1(:)).*kw2);b = std(phi(mask1(:)).*kw2);
imagesc(n1,m1,phi.*kw2.*mask1,[a-3*b a+3*b])
set(gca,'DataAspectRatio',[1 1 1]')
%set(gca,'YTick',0);set(gca,'XTick',0);set(gca,'TickLength',[0 0]);
t = colorbar('SouthOutside','peer',gca);
set(get(t,'ylabel'),'String', 'rad','FontName','Times','FontSize',14);
%scalebar('Location','southeast','Colour',barcolor,'Unit','um')
xlabel('(mm)'); ylabel('(mm)');
title(['Phase error - RMS = ' num2str(std(phi(mask1(:)).*kw2)) ' rad'] )


subplot(2,3,3)
a = mean(thickErr(mask1(:)));b = std(thickErr(mask1(:)));
imagesc(n1,m1,thickErr.*mask1,[a-3*b a+3*b])
set(gca,'DataAspectRatio',[1 1 1]')
%set(gca,'YTick',0);set(gca,'XTick',0);set(gca,'TickLength',[0 0]);
t = colorbar('SouthOutside','peer',gca);
set(get(t,'ylabel'),'String', 'um','FontName','Times','FontSize',14);
%scalebar('Location','southeast','Colour',barcolor,'Unit','um')
xlabel('(mm)'); ylabel('(mm)');
title(['Lens thickness error - RMS = ' num2str(std(thickErr(mask1(:)))) ' um'] )
% 
% 

subplot(2,3,4)
a = mean(phi1(mask1(:)));b = std(phi1(mask1(:)));
imagesc(n1,m1,phi1.*mask1.*1e9,[a-3*b a+3*b].*1e9)
set(gca,'DataAspectRatio',[1 1 1]')
%set(gca,'YTick',0);set(gca,'XTick',0);set(gca,'TickLength',[0 0]);
t = colorbar('SouthOutside','peer',gca);
set(get(t,'ylabel'),'String', 'nm','FontName','Times','FontSize',14);
%scalebar('Location','southeast','Colour',barcolor,'Unit','um')
xlabel('(mm)'); ylabel('(mm)');
title(['Wavefront error - RMS = ' num2str(std( phi1(mask1(:))*1e9 )) ' nm'] )


subplot(2,3,5)
a = mean(phi1(mask1(:)).*kw2);b = std(phi1(mask1(:)).*kw2);
imagesc(n1,m1,phi1.*kw2.*mask1,[a-3*b a+3*b])
set(gca,'DataAspectRatio',[1 1 1]')
%set(gca,'YTick',0);set(gca,'XTick',0);set(gca,'TickLength',[0 0]);
t = colorbar('SouthOutside','peer',gca);
set(get(t,'ylabel'),'String', 'rad','FontName','Times','FontSize',14);
%scalebar('Location','southeast','Colour',barcolor,'Unit','um')
xlabel('(mm)'); ylabel('(mm)');
title(['Phase error - RMS = ' num2str(std(phi1(mask1(:)).*kw2)) ' rad'] )


subplot(2,3,6)
a = mean(thickErr1(mask1(:)));b = std(thickErr1(mask1(:)));
imagesc(n1,m1,thickErr1.*mask1,[a-3*b a+3*b])
set(gca,'DataAspectRatio',[1 1 1]')
%set(gca,'YTick',0);set(gca,'XTick',0);set(gca,'TickLength',[0 0]);
t = colorbar('SouthOutside','peer',gca);
set(get(t,'ylabel'),'String', 'um','FontName','Times','FontSize',14);
%scalebar('Location','southeast','Colour',barcolor,'Unit','um')
xlabel('(mm)'); ylabel('(mm)');
title(['Lens thickness error - RMS = ' num2str(std(thickErr1(mask1(:)))) ' um'] )


% 
% %% ----------------------figure thickness error---------------------------%
% %% ----------------figure wavefront gradient with mask -------------------%
% figure(5)
% set(5,'units','points','position',[1550,1,650,400])
% 
% colormap(flipud(CM4))
% ma1 = mask1(1:3:end,1:3:end);td = thickErr(1:3:end,1:3:end);td(~ma1) = nan;
% [XZ,YZ] = meshgrid((1:1:size(td,2)).*pixsize./1000.*3,(1:1:size(td,1)).*pixsize./1000.*3);
% 
% 
% surf(XZ,YZ,td)
% shading interp
% axis tight
% view(-70,45)
% %camlight right
% xlabel('mm');ylabel('mm');zlabel('um');
% title(['Lens thickness error - RMS = ' num2str(std(thickErr(mask1(:)))) ' um'] )
% set(5,'PaperPositionMode','auto','PaperType','A3','PaperOrientation','landscape') 