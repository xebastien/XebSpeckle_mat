function [img3] = undistorImg(img,disto)


ROIdis = disto.ROI;
ROIdis(2) = ROIdis(1) + size(disto.h,1)-1;
ROIdis(4) = ROIdis(3) + size(disto.h,2)-1;
img1 = img(ROIdis(1):ROIdis(2),ROIdis(3):ROIdis(4));


[X,Y] = meshgrid(ROIdis(3):ROIdis(4),ROIdis(1):ROIdis(2));

U1 = disto.h + X;
V1 = disto.v + Y;

img2 = interp2(X,Y,single(img1), U1, V1,'cubic',mean(img1(:)));

img3 = single(img);
img3(ROIdis(1):ROIdis(2),ROIdis(3):ROIdis(4)) = img2;
