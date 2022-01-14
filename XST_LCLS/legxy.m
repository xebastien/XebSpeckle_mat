function [ modes, reconst] = legxy( sx,mask )
% Description: Represent a wavefront as a sum of Legendre polynomials using
%              a matrix inversion.
% 
% This function attempts to solve the a_i's in equation,
% 
%                     M
%                     __
%                    \
%       phi'      =  /__  a_i * P_i'(x,y)
%                    i=1
% 
% where the P_i'(x,y)'s are thederivative of the Legendre polynomials 
% phi is the wavefront to be represented as a sum of Legendre 
% polynomials, the a_i's are the legendre coefficients, and M is the number
% of legendrepolynomials to use. (up to 14)
%
% Input:    sx, sy - wavefront slopes to be represented as a sum of  legendre derivative polynomials
%                 that must be an nXmn array (rectangular)
%           (optional) M - Number of legendre polynomials to use (Default = 12)
% Output:   modes - legendre coefficients (a_i's) as a vector
%               reconst : wavefront reconstruction with the mode extracted    
% S.Berujon 


    M = 3;
if M > 19
    error('20? Why so high?.');
end



% catch working size
[nx, ny] = size(sx);

  
%% built a mask
  nn = max([nx ny]); 
if nargin == 2,
  
    mask = mask>0;    
else
   
  if nx == ny
      mask = ones(nn);
  elseif nn == nx
      mask = [ones(nx,ny) zeros(nx,nx-ny)];
      mask = circshift(mask,[0 floor((nx-ny)/2)]);
  elseif nn == ny
      mask = [ones(nx,ny) ; zeros(ny-nx,ny)];
      mask = circshift(mask,[floor((ny-nx)/2) 0]);
  end;  
end; 

%% load legendre polynom with normalization factor
F{1} = 'x' ;         n{1} = 3;      
F{2} = 'y ';         n{2} = 3;     
% F{3} = '3*x.^2-1' ;  n{3} = 5/4;    
% F{4} = '3*y.^2-1' ;  n{4} = 5/4;    
F{3} = 'x.*y' ;      n{3} = 9;     
% F{6} = '(3*x.^2-1).*y' ;    n{6} = 15/4;    dFx{6} = '6.*x.*y';        dFy{6} = '3*x.^2-1';
% F{7} = '(3*y.^2-1).*x' ;    n{7} = 15/4;    dFx{7} = '3*y.^2-1';       dFy{7} = '6*x.*y';
% F{8} = '(5*x.^2-3).*x' ;    n{8} = 7/4;     dFx{8} = '15.*x.^2-3';     dFy{8} = '0.*x';
% F{9} = '(5*y.^2-3).*y' ;    n{9} = 7/4;     dFx{9} = '0.*x';           dFy{9} = '15*y.^2-3';
% F{10} = '(35.*x.^4-30.*x.^2+3)' ;    n{10} = 9/64;   dFx{10} = '140.*x.^3-60.*x';                dFy{10} = '0.*x';
% F{11} = '(5*x.^2-3).*x.*y' ;         n{11} = 7/4*3;  dFx{11} = 'y.*(5.*x.^2-3)+10.*x.^2.*y';     dFy{11} = 'x.*(5.*x.^2-3)';
% F{12} = '(3*y.^2-1).*(3*x.^2-1)';    n{12} = 25/16;  dFx{12} = '6.*x.*(3.*y.^2-1)';              dFy{12} = '6.*y.*(3.*x.^2-1)';
% F{13} = '(5*y.^2-3).*y.*x' ;         n{13} = 7/4*3;  dFx{13} = 'y.*(5.*y.^2-3)';                 dFy{13} = 'x.*(5.*y.^2-3)+10.*y.^2.*x';
% F{14} = '(35.*y.^4-30.*y.^2+3)' ;    n{14} = 9/64;   dFx{14} = '0.*x';                           dFy{14} = '140.*y.^3-60.*y';
% 
% F{15} = '(63.*x.^5 - 70.*x.^3 + 15.*x)'; n{15} = 11/64;     dFx{15} = '315*x.^4 - 210*x.^2 + 15';       dFy{15} = '0.*y';
% F{16} = '(63.*y.^5 - 70.*y.^3 + 15.*y)'; n{16} = 11/64;     dFx{16} = '0.*x';       dFy{16} = '315*y.^4 - 210*y.^2 + 15';
% 
% F{17} = '(231.*x.^6 - 315.*x.^4 + 105.*x.^2-5)'; n{17} = 13/256;     dFx{17} = '1386*x.^5 - 1260*x.^3 + 210.*x';       dFy{17} = '0.*y';
% F{18} = '(231.*y.^6 - 315.*y.^4 + 105.*y.^2-5)'; n{18} = 11/256;     dFx{18} = '0.*x';       dFy{18} = '1386*y.^5 - 1260*y.^3 + 210.*y';
% 
% 
% F{19} = '1' ;        n{19} = 1;     dFx{19} = '0.*x';    dFy{19} = '0.*x';



%% create a normalized size
X = -1:1/(nn/2-1/2):1;
Y = -1:1/(nn/2-1/2):1;

[x2,y2] = meshgrid(X,Y);

x = x2(mask == 1);
y = y2(mask == 1);

%% calculate the modes values
B = zeros(numel(x),M);


for kk = 1:M
    B(:,kk) = reshape(eval(F{kk}).*sqrt(n{kk}),[],1) ;
end;

%modes = pinv(B)*([reshape(sx(mask == 1),[],1) ; reshape(sy(mask == 1),[],1)]);
modes = pinv(B)*(reshape(sx(mask == 1),[],1));
modes = modes' ;%.* nn/2;



%% reconstruct with the tilt included
reconst = zeros(nn,nn);
x=x2;y=y2;
for kk = 3: M 
    reconst = reshape(eval(F{kk}).*sqrt(n{kk}).*modes(kk),nn,nn) + reconst; 
end;
reconst = reconst - mean(reconst(:));

reconst = reconst(1:nx,1:ny);
% %% reconstruct without the tilt
% 
% reconst_wt = zeros(nx,ny);
% 
% for kk = 3:M
%     reconst_wt = reshape(eval(F{kk}).*sqrt(n{kk}).*modes(kk+1),nx,ny) + reconst_wt; 
% end;























