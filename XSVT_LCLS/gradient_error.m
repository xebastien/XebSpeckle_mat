function [reconst1, reconst2] = gradient_error(SpeckDisH,SpeckDisV,mask)
SpeckDisH = SpeckDisH-mean(SpeckDisH(:));
SpeckDisV = SpeckDisV-mean(SpeckDisV(:));
% % pixsize in meters
% 
% 
% [m n] = size(SpeckDisH);
% 
% % horizontal
% x = (1:1:m) - round(mean(m/2));
% leg = repmat(x,n,1);   
% [m n] = size(SpeckDisH);    
% phsp = reshape(SpeckDisH,m*n,1);
% % (3^(1/2)).*
% Z = reshape(leg,m*n,1);
% a_h = pinv(Z)*phsp;
% reconst1 = a_h *leg;
% 
% [m n] = size(SpeckDisV);
% % horizontal
% x = (1:1:m)' - round(mean(m/2));
% x = x.*5.8e-6;
% leg = repmat(x,1,n);   
% [m n] = size(SpeckDisV);    
% phsp = reshape(SpeckDisV,m*n,1);
% % (3^(1/2)).*
% Z = reshape(leg,m*n,1);
% a_h = pinv(Z)*phsp;
% reconst2 = a_h *reshape(leg,m,n);

% return;

if nargin == 2,
    mask = ones(size(SpeckDisH));
end;
mask = mask>0;
[m, n] = size(SpeckDisH);    

x = linspace(-1,1,n);
leg = repmat(x,m,1);   
phsp = SpeckDisH(mask(:) == 1);
a_h = pinv(leg(mask(:)==1))*phsp;
reconst1 = a_h *leg;

% vertical
% x = -1:1/(128-1/2):1;
% leg= repmat(x,256,1);   
% [m n] = size(rot90(SpeckDisV,1));    
% phsp = reshape(rot90(SpeckDisV,1),m*n,1);
% 
% Z = reshape(imresize(leg,[m n]),m*n,1);
% a_h = pinv(Z)*phsp;
% reconst2 = rot90(a_h *imresize(leg,[m n]),-1);



x = linspace(-1,1,m);
leg = repmat(x,n,1);  
SpeckDisV = rot90(SpeckDisV,1);

mask = rot90(mask,1);
phsp = SpeckDisV(mask(:));
a_h = leg(mask(:))\phsp;
reconst2 = a_h *leg;

reconst2 = rot90(reconst2,-1);
% SpeckDisV = rot90(SpeckDisV,-1);
%  mask=rot90(mask,-1);
return;
