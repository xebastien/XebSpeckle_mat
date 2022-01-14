function [delay,coeff,pic3] = SliceDelayAbsMode(plane1, plane2)

%calculate the delay in step between two slices made from the third
%dimension of a stack images.

[m, n] = size(plane1);
%pic = [0 0 0 0 0];

optproc = 1;% 1 for my linux 2- for windows people
optdisp = 1;% to activate display of the peak



edge = 50 ;%important for the offset
edgeV = 500;%
plane1 = plane1./repmat(sum(plane1,2),[1 n])./repmat(sum(plane1,1),[m 1]);
plane2 = plane2./repmat(sum(plane2,2),[1 n])./repmat(sum(plane2,1),[m 1]);

plane1 = (plane1 - mean(plane1(:)))./std(plane1(:));
plane2 = (plane2 - mean(plane2(:)))./std(plane2(:));
plane1t = plane1(edgeV:end-edgeV,edge:end-edge);


if optproc == 1
    [crosscorr3,pic3] = normxcorr2_mexsub(plane1t,plane2,'valid');
    xpeak = pic3(3); ypeak = pic3(4);
    delay = [xpeak ypeak];
    coeff  = pic3(5);
    crosscorr = crosscorr3;
end;


if optproc == 2
   crosscorr3 = normxcorr2(plane1t,plane2);
   [xpeak, ypeak, coeff] = findpeak(crosscorr3,true);
   corr_offset = [ (ypeak-size(plane1t,1)-edgeV+1)  (xpeak-size(plane1t,2)-edge+1)];
   delay = corr_offset;
   pic3 = [corr_offset(1) corr_offset(2) corr_offset(1) corr_offset(2) coeff];
   crosscorr = crosscorr3(m-2*edgeV:m,n-2*edge:n);
end;
flag = 0;
if (pic3(4) > -1.5 )&& (pic3(4) < 1.5)
    gaussEdge =3;
    % be carefull of which max2 function is used(id19 function renders
    % three values)
    [~,I(1),I(2)] = max2(crosscorr);
    yd = crosscorr(I(1),I(2)-3:I(2)+3);
    x1d = I(2)-gaussEdge:1:I(2)-gaussEdge+length(yd)-1;
    sigma0 = .3;

    f = @(xl) sum(((xl(4).*exp( -(x1d - xl(1)).^2 / (2*xl(2)^2))+xl(3))-yd).^2);
    inp = [x1d(gaussEdge+1) sigma0 0 max(yd)];
    out = fminsearch(f,inp);
    
    if abs(out(1)) < (edge-3), delay(2) = out(1)-edge;  end;
    y2 = out(4).*exp( -(x1d - out(1)).^2 / (2*out(2)^2))+out(3);
    flag =1;

%     [cx,cy,sx,sy,PeakOD] = Gaussian2Dfit(crosscorr(:,I(2)-gaussEdge:I(2)+gaussEdge),1e-6);
%     cx2 = cx - gaussEdge;
%     [sizey, sizex] = size(crosscorr(:,I(2)-6:I(2)+6));
%     [x,y] = meshgrid(1:sizex,1:sizey);
%     fitout = abs(PeakOD)*(exp(-0.5*(x-cx).^2./(sx^2)-0.5*(y-cy).^2./(sy^2)));

end;




if optdisp
    %figure(1)
    subplot(2,2,1)
    imagesc(plane1t)

    subplot(2,2,2)
    imagesc(plane2)

    subplot(2,2,3)
    %surf(crosscorr)
    %shading interp
    imagesc(crosscorr)

    subplot(2,2,4)
    [~,I] = max2(crosscorr);
    plot(crosscorr(I(1),:))
    if flag
    hold on
    plot(x1d,y2,'r')
    hold off



%     figure(2)
%     surf(crosscorr(:,I(2)-gaussEdge:I(2)+gaussEdge))
%     shading interp
%     alpha(.5)
%     hold on
%     surf(fitout)
%     shading interp
%     hold off

    end;
    drawnow
end;






