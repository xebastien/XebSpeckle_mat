CM1 = colormap('gray');% colormap for global phsae graident
CM2 = colormap(othercolor('PuBu9'));%cbrewer('seq','PuBu',100);%;colormap(othercolor('Reds9'));%cbrewer('div','RdBu',100);%cbrewer('div','RdGy',100);%cbrewer('div','PuOr',100);%
%colormap(othercolor('RdYlBu11'));% colormap for global phase gradient error
CM4 = cbrewer('div','RdGy',100);%colormap(othercolor('PRGn5'));% colormap for phase error and thickness error
CM3 = cbrewer('div','RdBu',100);%colormap('gray');
barcolor = [1 0 0];

%% ----------------figure wavefront gradient as calculated---------------%
dphaseX = vert2;dphaseY = horiz2;
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
nbstd = 2;
subplot(1,2,1)
a = mean(dphaseX(:).*1e6);b = std(dphaseX(:).*1e6);
imagesc(n1,m1,dphaseX.*1e6,[a-nbstd*b a+nbstd*b])
set(gca,'DataAspectRatio',[1 1 1]')
%set(gca,'YTick',0);set(gca,'XTick',0);set(gca,'TickLength',[0 0]);
t = colorbar('SouthOutside','peer',gca);
set(get(t,'ylabel'),'String', 'urad','FontName','Times','FontSize',14);
scalebar('Location','southeast','Colour',barcolor,'Unit','mm','ScaleLength',round(m1(end)/2*10)/1e1)
xlabel(['(mm)  - Lens center at ' num2str(fliplr(Center2(2)))]);
ylabel(['(mm)  - Lens center at ' num2str(fliplr(Center2(1)))]);
title('Vertical wavefront gradient')
grid off


subplot(1,2,2)
a = mean(dphaseY(:).*1e6);b = std(dphaseY(:).*1e6);
imagesc(n1,m1,dphaseY.*1e6,[a-nbstd*b a+nbstd*b])
set(gca,'DataAspectRatio',[1 1 1]')
%set(gca,'YTick',0);set(gca,'XTick',0);set(gca,'TickLength',[0 0]);
t = colorbar('SouthOutside','peer',gca);
set(get(t,'ylabel'),'String', 'urad','FontName','Times','FontSize',14);
scalebar('Location','southeast','Colour',barcolor,'Unit','um','ScaleLength',round(m1(end)/2*10)/1e1)
xlabel(['(mm)  - Lens center at ' num2str(fliplr(Center2(2)))]);
ylabel(['(mm)  - Lens center at ' num2str(fliplr(Center2(1)))]);
title('Horizontal wavefront gradient')
grid off

set(1,'PaperPositionMode','auto','PaperType','A3','PaperOrientation','landscape')  


%% ----------------figure wavefront gradient with mask -------------------%
[m,n] = size(xg);
m1 = (1:1:m).*pixsize./1000;n1 = (1:1:n).*pixsize./1000;
m1 = m1 - mean(m1);n1 = n1 - mean(n1);

figure(2)
colormap(flipud(CM2))
set(2,'units','points','position',[2000,400,800,400])
nbstd = 3;

subplot(1,2,1)
a = 0;b = std(xg(mask1(:)));
imagesc(n1,m1,xg.*mask1 - mean(xg(mask1(:))),[a-nbstd*b a+nbstd*b])
set(gca,'DataAspectRatio',[1 1 1]')
%set(gca,'YTick',0);set(gca,'XTick',0);set(gca,'TickLength',[0 0]);
t = colorbar('SouthOutside','peer',gca);
set(get(t,'ylabel'),'String', 'urad','FontName','Times','FontSize',14);
%scalebar('Location','southeast','Colour',barcolor,'Unit','um')
xlabel(['(mm) - Focal length = ' num2str(raddiff1) ' m']);
ylabel('(mm)');
title('Horizontal wavefront gradient error')

subplot(1,2,2)
a = 0;b = std(yg(mask1(:)));
imagesc(n1,m1,yg.*mask1- mean(yg(mask1(:))),[a-nbstd*b a+nbstd*b])
set(gca,'DataAspectRatio',[1 1 1]')
%set(gca,'YTick',0);set(gca,'XTick',0);set(gca,'TickLength',[0 0]);
t = colorbar('SouthOutside','peer',gca);
set(get(t,'ylabel'),'String', 'urad','FontName','Times','FontSize',14);
%scalebar('Location','southeast','Colour',barcolor,'Unit','um')
xlabel(['(mm) - Focal length = ' num2str(raddiff2) ' m']);
ylabel('(mm)');
title('Vertical wavefront gradient error')

set(2,'PaperPositionMode','auto','PaperType','A3','PaperOrientation','landscape')  
%% ---------------------------figure masks -------------------------------%
figure(3)
set(3,'units','points','position',[1800,200,800,400])

subplot(1,3,1)
b1 = 3*std(reshape(vert2,1,[]));
imagesc((~mask1).*vert2,[-b1 b1])
set(gca,'DataAspectRatio',[1 1 1]')
set(gca,'YTick',0);set(gca,'XTick',0);set(gca,'TickLength',[0 0]);
title(['Mask diameter = ' num2str(maskDiam)]);
xlabel(['Center location in pixel ' num2str(fliplr(Center))]);

subplot(1,3,2)
imagesc(file1./max(file1(:))+mask2.*0.5)
set(gca,'DataAspectRatio',[1 1 1]')
set(gca,'YTick',0);set(gca,'XTick',0);set(gca,'TickLength',[0 0]);
title(['Mask diameter = ' num2str(maskDiam)]);
xlabel(['Center location in pixel ' num2str(fliplr(Center))]);

subplot(1,3,3)
imagesc(file1)
set(gca,'DataAspectRatio',[1 1 1]')
set(gca,'YTick',0);set(gca,'XTick',0);set(gca,'TickLength',[0 0]);
title(['Mask diameter = ' num2str(maskDiam)]);
xlabel(['Center location in pixel ' num2str(fliplr(Center))]);

colormap(CM3)

set(3,'PaperPositionMode','auto','PaperType','A3','PaperOrientation','landscape')  


%% --------------figure wavefront error and thickness error --------------%
%% ----------------figure wavefront gradient with mask -------------------%
[m,n] = size(phi);
m1 = (1:1:m).*pixsize./1000;n1 = (1:1:n).*pixsize./1000;
m1 = m1 - mean(m1);n1 = n1 - mean(n1);

figure(4)
set(4,'units','points','position',[ 420.0000 ,173.2500, 942.0000, 576.7500])

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
xlabel(['(mm) PV = ' num2str(min(phi1(:).*kw2.*mask1(:))) '  ' num2str(max(phi1(:).*kw2.*mask1(:)))]); ylabel('(mm)');
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
title(['T. error , no  defoc xy- RMS = ' num2str(std(thickErr1(mask1(:)))) ' um'] )

set(4,'PaperPositionMode','auto','PaperType','A3','PaperOrientation','landscape') 
return;
%% export drawings
writeFolder = '/data/bm05/imaging/seb/LCLS/processing/XSVT/plate';
writeFolder = '/data/bm05/imaging/seb/LCLS/processing/XSVT/noplate';

print(2,'-dpdf',fullfile(writeFolder,'Wavefront_gradient_error.pdf'));
saveas(2, fullfile(writeFolder,'Wavefront_gradient_error.svg'));
print(2,'-dpng',fullfile(writeFolder,'Wavefront_gradient_error.png'));

print(4,'-dpdf',fullfile(writeFolder,'Wavefront_error.pdf'));
saveas(4, fullfile(writeFolder,'Wavefront_error.svg'));
print(4,'-dpng',fullfile(writeFolder,'Wavefront_error.png'));

save(fullfile(writeFolder,'phi'),'phi1')

return;
%plot cut to compare - need both phi1 in memory
mil = round(size(phi1)/2);
figure(5)
subplot(2,1,1)
plot(n1,phiplate(mil,:).*kw2.*mask1(mil,:),'b')
hold on

subplot(2,1,1)
plot(n1,phinoplate(mil,:).*kw2.*mask1(mil,:),'r')
xlabel('mm');ylabel('radian')
legend('With plate', 'Without plate')

subplot(2,1,2)
plot(m1,phiplate(:,mil).*kw2.*mask1(:,mil),'b')
hold on

subplot(2,1,2)
plot(m1,phinoplate(:,mil).*kw2.*mask1(:,mil),'r')
xlabel('mm');ylabel('radian')
legend('With plate', 'Without plate')





return






%% ----------------------figure thickness error---------------------------%
%% ----------------figure wavefront gradient with mask -------------------%
figure(5)
set(5,'units','points','position',[1550,1,650,400])

colormap(flipud(CM3))
ma1 = mask1(1:5:end,1:5:end);
%td = thickErr(1:9:end,1:9:end);
td = -phiFull(1:15:end,1:15:end);
%td(~ma1) = nan;
[XZ,YZ] = meshgrid((1:1:size(td,2)).*pixsize./1000.*3,(1:1:size(td,1)).*pixsize./1000.*3);


surf(XZ,YZ,td.*1e9)
%shading interp
axis tight
view(-70,45)
%camlight right
xlabel('mm');ylabel('mm');zlabel('nm');
title(['Lens thickness error - RMS = ' num2str(std(thickErr(mask1(:)))) ' um'] )
set(5,'PaperPositionMode','auto','PaperType','A3','PaperOrientation','landscape') 