%  SolveLSChol.m
% frioom the net
%  Use sparse least squares reconstructor to estimate pupil-plane phase.
function wft = WftSolveLSChol(sx_d,sy_d,h,MASK)
  
sx_d = sx_d.*2*h;  
sy_d = sy_d.*2*h;

if nargin < 4 , MASK = ones(size(sx_d));end;

[nx, ny] = size(sx_d);
  
nn = max([nx ny]);
  
if nx == ny
      
    mask = ones(nn);
      
      
elseif nn == nx
    sx_d = [sx_d zeros(nx,nx-ny)];
    sy_d = [sy_d zeros(nx,nx-ny)];
      
    mask = [ones(nx,ny) zeros(nx,nx-ny)];
    mask = [MASK zeros(nx,nx-ny)];
      
elseif nn == ny
    sx_d = [sx_d ; zeros(ny-nx,ny)];
    sy_d = [sy_d ; zeros(ny-nx,ny)];
      
      
    mask = [ones(nx,ny) ; zeros(ny-nx,ny)];
    mask = [MASK ; zeros(ny-nx,ny)];  
end;  
  
M = spdiags(mask(:), 0, nn^2,nn^2);
  
[Gx,Gy,Gamma_x_hat,Gamma_y_hat] = mk_GammaFried2D2(nn,nn);

  
% Enter the coefficient matrix
A = Gx'*M*Gx + Gy'*M*Gy+sqrt(eps)*speye((nn+1)^2,(nn+1)^2);
% Compute its Cholesky factorization with sparse reordering.
I = speye(size(A));
p = symamd(A);
C = chol(A(p,p));
Ip = I(:,p);                   % Inverse permutation matrix.
% Right-hand side vector.
b = Gx'*M*sx_d(:) + Gy'*M*sy_d(:);
% Solution via Cholesky Factorization.
phi_LS = Ip* ( C \ ( C' \ b(p) ) );
phi_LS = reshape(phi_LS,nn+1,nn+1);


%   figure(2)
%      subplot(221)
%       Mphi_LS = phi_LS(1:nx,1:ny);
%       imagesc(Mphi_LS), colorbar, colormap(hot)
%       title('LS Reconstruction of Phase')
%    
  wft = phi_LS(1:nx,1:ny);
  
% imagesc(Gx(end-1000:end,end-1000:end))