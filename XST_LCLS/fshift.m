function [x_offset,y_offset] = fshift(a,t)


% Determine padding size in x and y dimension
size_t      = size(t);
size_a      = size(a);
outsize     = size_t + size_a - 1;

% Determine 2D cross correlation in Fourier domain
Ft = fft2(t, outsize(1), outsize(2));
Fa = fft2(a, outsize(1), outsize(2));
c = abs( fftshift( ifft2(Fa .* conj(Ft))) );

% Find peak
[max_c, imax]   = max(abs(c(:)));
[ypeak, xpeak]  = ind2sub(size(c), imax(1));

% Correct found peak location for image size
corr_offset = round([(ypeak-(size(c, 1)+1)/2) (xpeak-(size(c, 2)+1)/2)]);

% Write out offsets
y_offset = corr_offset(1);
x_offset = corr_offset(2);