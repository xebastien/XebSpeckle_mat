function [Gamma_x,Gamma_y,Gamma_x_hat,Gamma_y_hat] = mk_GammaFried2D2(nx,ny)
%  [Gamma_x,Gamma_y,Gamma_x_hat,Gamma_y_hat] = mk_GammaFried(nx)
%
%  Set up spatial domain sparse matrix representors for 
%  x- and y-components of discrete gradient operator corresponding 
%  to Fried geometry. Then set up Fourier domain representors.


%  Compute spatial domain representors.
  onevecx = ones(nx+1,1);
  onevecy = ones(ny+1,1);
  G1 = spdiags([-onevecx onevecx], [0 1], nx,nx+1);
  G2 = spdiags([onevecy onevecy], [0 1], ny+1,ny+1);
% G1(nx,1) = 1;
%  G2(ny,1) = 1;
  Gamma_x = kron(G2,G1);
  Gamma_x = Gamma_x(1:(nx*ny),:);
  G1 = spdiags([-onevecy onevecy], [0 1], ny+1,ny+1);
  G2 = spdiags([onevecx onevecx], [0 1], nx,nx+1);
 %G1(ny,1) = 1;
 % G2(nx,1) = 1;
  Gamma_y = kron(G1,G2);
  Gamma_y = Gamma_y(1:(nx*ny),:);
% %  Set up Fourier domain representors. 
  Sx_hat = fourier_xshift(nx,ny);  %  Shift in x-direction
  Sy_hat = fourier_yshift(nx,ny);  %  Shift in y-direction
  Gamma_x_hat = (Sx_hat - 1).*(Sy_hat + 1);
  Gamma_x_hat = spdiags(Gamma_x_hat(:), 0, nx*ny,nx*ny);
  Gamma_y_hat = (Sx_hat + 1).*(Sy_hat - 1);
  Gamma_y_hat = spdiags(Gamma_y_hat(:), 0, nx*ny,nx*ny);
%figure(5);subplot(3,1,1);imagesc(G1);subplot(3,1,2);imagesc(G2);subplot(3,1,3);imagesc(Gamma_y);


 function Sx_hat = fourier_xshift(nx,ny)

%  Sx_hat = fourier_xshift(nx,ny)
%
%  Compute Fourier representer Sx_hat for operator Sx for which
%     [Sx*F](i,j) = F(i+1,j),
%  where F is periodically extended nx X ny array. 
%  To compute Sx*F, evaluate 
%     real(ifft2(Sx_hat.*fft2(F))).

  kappa_x = (0:nx-1)'/nx;
  imath = sqrt(-1);
  kappa = kappa_x*ones(1,ny);
  Sx_hat = exp(imath*2*pi*kappa);
 
  
   function Sy_hat = fourier_yshift(nx,ny)

%  Sy_hat = fourier_yshift(nx,ny)
%
%  Compute Fourier representer Sy_hat of operator Sy for which
%     [Sy*F](i,j) = F(i,j+1),
%  where F is periodically extended nx X ny array. 
%  To compute Sy*F, evaluate 
%     real(ifft2(Sy_hat.*fft2(F))).

  kappa_y = (0:ny-1)/ny;
  imath = sqrt(-1);
  kappa = ones(nx,1)*kappa_y;
  Sy_hat = exp(imath*2*pi*kappa);
