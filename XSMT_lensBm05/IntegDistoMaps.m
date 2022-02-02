
%pixsizeH = pixsizeV;
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      Integrate the 2x2 derivative disto maps to get disto maps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



sdhh1 = sdhh;
sdhh1(abs(sdhh1 - mean(sdhh1(maskI))) >(5*std(sdhh1(:)))) = mean(sdhh1(maskI));
sdvh1 = sdvh;
sdvh1(abs(sdvh1 - mean(sdvh1(maskI))) >(5*std(sdvh1(:)))) = mean(sdvh1(maskI));
sdvv1 = sdvv;
sdvv1(abs(sdvv1 - mean(sdvv1(maskI))) >(5*std(sdvv1(:)))) = mean(sdvv1(maskI));
sdhv1 = sdhv;
sdhv1(abs(sdhv1 - mean(sdhv1(maskI))) >(5*std(sdhv1(:)))) = mean(sdhv1(maskI));

sdhh1 = sdhh1 - mean(sdhh1(maskI));sdvh1 = sdvh1 - mean(sdvh1(maskI));
sdvv1 = sdvv1 - mean(sdvv1(maskI));sdhv1 = sdhv1 - mean(sdhv1(maskI));

sdhh1(~maskI) = 0;sdvh1(~maskI) = 0;
sdvv1(~maskI) = 0;sdhv1(~maskI) = 0;

disth1 = intgrad2(sdhh1,sdvh1);
distv1 = intgrad2(sdhv1,sdvv1);

disth2 = frankotchellappa(sdhh1,sdvh1);
distv2 = frankotchellappa(sdhv1,sdvv1);

disth3 = WftSolveLSChol(sdvh1,sdhh1,1,ones(size(sdhh1)));
distv3 = WftSolveLSChol(sdvv1,sdhv1,1,ones(size(sdhh1)));

disth = disth1 - mean(disth1(maskI));
distv = distv1 - mean(distv1(maskI));


disth1 = disth1 - mean(disth1(maskI));
distv1 = distv1 - mean(distv1(maskI));

disth2 = disth2 - mean(disth2(maskI));
distv2 = distv2 - mean(distv2(maskI));

disth3 = disth3 - mean(disth3(maskI));
distv3 = distv3 - mean(distv3(maskI));
% save('','disth','distv')

a = min([sdhh(:) ;sdvv(:); sdhv(:) ;sdvh(:)]);
b = max([sdhh(:); sdvv(:) ;sdhv(:) ;sdvh(:)]);
figure(4)
subplot(2,2,1)
imagesc(sdhh)
title('Dx delta u ')
xlabel('pixels')
ylabel('pixels')
colorbar
subplot(2,2,2)
imagesc(sdhv)
xlabel('pixels')
ylabel('pixels')
title('Dx delta v ')
colorbar
subplot(2,2,3)
imagesc(sdvh)
xlabel('pixels')
ylabel('pixels')
title('Dy delta u')
colorbar
subplot(2,2,4)
imagesc(sdvv)
xlabel('pixels')
ylabel('pixels')
title('Dy delta v')
colorbar
colormap(othercolor('Paired4'))


figure(2)
subplot(1,2,1)
imagesc(disth)
title('Delta u')
xlabel('pixels')
ylabel('pixels')
colorbar
subplot(1,2,2)
imagesc(distv)
title('Delta v')
xlabel('pixels')
ylabel('pixels')
colorbar
colormap(othercolor('RdYlGn10'))
colormap(othercolor('Paired4'))



disto.v = distv;
disto.h = disth;
%return;
save(saveplace,'disto');
save(saveplace,'ROI','pixsizeV','pixsizeH','-append');
   
