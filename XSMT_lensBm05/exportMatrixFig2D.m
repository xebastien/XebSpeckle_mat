legendshere = cell(1,length(subFolders));
for k1 = 1:length(subFolders),legendshere{k1} = strrep(strrep(subFolders(k1).name, '.mat',''),'_','-');end;


nbstd = 2;

figure(5)
set(5,'units','points','position',[1600,450,1200,600])


subplot(2,2,1)
p = plot(diamet(1:end), focV(:,1:end));
if size(focM,1)>5,set(p(6:end),'LineStyle','-.');end;
a = mean(focV(:));b = std(focV(:));
axis([0 diamet(1) (a-nbstd *b) (a+nbstd *b)])
legend(legendshere {:},'Location','EastOutside')%,'Orientation','horizontal') 
xlabel('Mask diameter (um)')
ylabel('Focal distance (m)')
title('Vertical focal length / mask diameter' );

subplot(2,2,2)
p = plot(diamet(1:end), focH(:,1:end));
if size(focM,1)>5,set(p(6:end),'LineStyle','-.');end;
a = mean(focH(:));b = std(focH(:));
axis([0 diamet(1) (a-nbstd *b) (a+nbstd *b)])
legend(legendshere {:},'Location','EastOutside')%,'Orientation','horizontal') 
xlabel('Mask diameter (um)')
ylabel('Focal distance (m)')
title('Horizontal focal length / mask diameter' );

subplot(2,2,3)
p = plot(diamet(1:end), focM(:,1:end).*2.*delta.*1e6);
if size(focM,1)>5,set(p(6:end),'LineStyle','-.');end;
a = mean(focM(:).*2.*delta.*1e6);b = std(focM(:).*2.*delta.*1e6);
axis([0 diamet(1) (a-nbstd *b) (a+nbstd *b)])
legend(legendshere {:},'Location','EastOutside')%,'Orientation','horizontal') 
xlabel('Mask diameter (um)')
ylabel('Radius (um)')
title('Effective radius / mask diameter' );


subplot(2,2,4)
p = plot(diamet(1:end), ter(:,1:end));
if size(focM,1)>5,set(p(6:end),'LineStyle','-.');end;
a = mean(ter(:));b = std(ter(:));
axis([0 diamet(1) (a-3*b) (a+3 *b)])
legend(legendshere {:},'Location','EastOutside')%,'Orientation','horizontal') 
xlabel('Mask diameter (um)')
ylabel('RMS thickness error  (um)')
title('RMS thickness error / mask diameter' );


set(5,'PaperPositionMode','auto','PaperType','A3','PaperOrientation','landscape')  