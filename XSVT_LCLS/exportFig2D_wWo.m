CM1 = colormap('gray');% colormap for global phsae graident
CM2 = colormap(othercolor('PuBu9'));%cbrewer('seq','PuBu',100);%;colormap(othercolor('Reds9'));%cbrewer('div','RdBu',100);%cbrewer('div','RdGy',100);%cbrewer('div','PuOr',100);%
%colormap(othercolor('RdYlBu11'));% colormap for global phase gradient error
CM4 = cbrewer('div','RdGy',100);%colormap(othercolor('PRGn5'));% colormap for phase error and thickness error
CM3 = cbrewer('div','RdBu',100);%colormap('gray');
barcolor = [1 0 0];


bmin = -3*pi;bmax = -bmin;
%% ----------------figure wavefront gradient as calculated---------------%
dphaseX = vert2;dphaseY = horiz2;
[m,n] = size(dphaseX);
m1 = (1:1:m).*pixsize./1000;n1 = (1:1:n).*pixsize./1000;

Center2 = Center.*pixsize./1000;
Center2(2) = Center2(2)- mean(n1);
Center2(1) = Center2(1)- mean(m1);
Center2 =  round(Center2.*1000)./1000;
m1 = m1 - mean(m1);n1 = n1 - mean(n1);







%% ------------------------cuts plots------------------_%
%plot cut to compare - need both phi1 in memory
phi1 = phiWo;

mil = round(size(phi1)/2);
phi1 = phi1 - mean(phi1(:));

%close(1)
figure(6)
%set(6,'units','points','position',[ 220.0000 ,273.2500, 942.0000, 326.7500])
subplot(2,3,2)
hold off
a = mean(phi1(mask1(:)).*kw2);b = std(phi1(mask1(:)).*kw2)./1.15;
imagesc(n1,m1,phi1.*mask1.*kw2,[bmin bmax])
hold on
xlabel('(mm)')

line((n1),(n1).*0,'color','black','LineStyle',':','LineWidth',2)
hold on
line((m1).*0,(m1),'color','black','LineStyle',':','LineWidth',2)
hold off

colormap((CM4))
grid on
set(gca,'DataAspectRatio',[1 1 1]')
%set(gca,'YTick',0);set(gca,'XTick',0);set(gca,'TickLength',[0 0]);
t = colorbar('WestOutside','peer',gca);
set(get(t,'ylabel'),'String', '(rad)','FontName','Times','FontSize',14);
%scalebar('Location','southeast','Colour',barcolor,'Unit','mm','ScaleLength',0.30)

%set(get(t,'ylabel'),'String', 'nm','FontName','Times','FontSize',14);
%mil = [0 0];
%grid on
hold off


%% ------------------------cuts plots------------------_%
%plot cut to compare - need both phi1 in memory
phi1 = phiW;

mil = round(size(phi1)/2);
phi1 = phi1 - mean(phi1(:));

%close(1)
figure(6)
%set(6,'units','points','position',[ 220.0000 ,273.2500, 942.0000, 326.7500])
subplot(2,3,1)
hold off
% a = mean(phi1(mask1(:)).*kw2);b = std(phi1(mask1(:)).*kw2);
imagesc(n1,m1,phi1.*mask1.*kw2,[bmin bmax])
hold on
xlabel('(mm)')

line((n1),(n1).*0,'color','b','LineStyle',':','LineWidth',2)
hold on
line((m1).*0,(m1),'color','b','LineStyle',':','LineWidth',2)
hold off


colormap((CM4))
grid on
set(gca,'DataAspectRatio',[1 1 1]')
%set(gca,'YTick',0);set(gca,'XTick',0);set(gca,'TickLength',[0 0]);
t = colorbar('WestOutside','peer',gca);
set(get(t,'ylabel'),'String', '(rad)','FontName','Times','FontSize',14);
%scalebar('Location','southeast','Colour',barcolor,'Unit','mm','ScaleLength',0.30)

%set(get(t,'ylabel'),'String', 'nm','FontName','Times','FontSize',14);
%mil = [0 0];
%grid on
hold off


%% ------------------------cuts plots------------------_%
%plot cut to compare - need both phi1 in memory
phi1 = phiD;

mil = round(size(phi1)/2);
phi1 = phi1 - mean(phi1(:));
[m2,n2] = size(phi1);

%close(1)
figure(6)
%set(6,'units','points','position',[ 220.0000 ,273.2500, 942.0000, 326.7500])
subplot(2,3,3)
hold off
% a = mean(phi1(mask1(:)).*kw2);b = std(phi1(mask1(:)).*kw2);
imagesc(n1,m1,phi1.*mask1.*kw2,[bmin bmax])
hold on
xlabel('(mm)')

line((n1),(n1).*0,'color','m','LineStyle',':','LineWidth',2)
hold on
line((m1).*0,(m1),'color','m','LineStyle',':','LineWidth',2)
hold off

colormap((CM4))
grid on
set(gca,'DataAspectRatio',[1 1 1]')
%set(gca,'YTick',0);set(gca,'XTick',0);set(gca,'TickLength',[0 0]);
t = colorbar('WestOutside','peer',gca);
set(get(t,'ylabel'),'String', '(rad)','FontName','Times','FontSize',14);
%scalebar('Location','southeast','Colour',barcolor,'Unit','mm','ScaleLength',0.30)

%set(get(t,'ylabel'),'String', 'nm','FontName','Times','FontSize',14);
%mil = [0 0];
%grid on
hold off






subplot(3,2,5)
% yyaxis left
% plot(n1,phi1(mil(1),:).*1e9.*mask1(mil(1),:),'r')
% xlabel('mm');ylabel('nm')
% yyaxis right
plot(n1,phiW(mil(1),:).*kw2.*mask1(mil(1),:),'b','LineWidth',2)
hold on
plot(n1,phiWo(mil(1),:).*kw2.*mask1(mil(1),:),'black','LineWidth',2)
plot(n1,phiD(mil(1),:).*kw2.*mask1(mil(1),:),'m','LineWidth',1,'LineStyle','--')
xlabel('Horizontal axis (mm)');ylabel('(rad)')
axis([-0.6 0.6 -12 6.3])
hold off

subplot(3,2,6)
% yyaxis left
% plot(m1,phi1(:,mil(2)).*1e9.*mask1(:,mil(2)),'g')
% xlabel('mm');ylabel('nm')
% yyaxis right
plot(n1,phiW(:,mil(2)).*kw2.*mask1(:,mil(2)),'b','LineWidth',2)%,'LineStyle',':')
hold on
plot(n1,phiWo(:,mil(2)).*kw2.*mask1(:,mil(2)),'black','LineWidth',2)%,'LineStyle','--')
plot(n1,phiD(:,mil(2)).*kw2.*mask1(:,mil(2)),'m','LineWidth',1,'LineStyle','--')
xlabel('Vertical axis (mm)');ylabel('(rad)')
%legend()
hold off
legend('With Plate','Without Plate','Differential','Orientation','horizontal')
axis([-0.6 0.6 -12 6.3])

leng = min(m,n);
cutsout = cat(1,m1(1:leng),phi((1:leng),mil(2))'.*1e9.*mask1((1:leng),mil(2))',phi(mil(1),(1:leng)).*1e9.*mask1(mil(1),(1:leng)))';
set(6,'PaperPositionMode','auto','PaperType','A3','PaperOrientation','landscape') 



