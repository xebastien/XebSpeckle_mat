function [a, reconst]= legendre_coeffs(phi,M)

% Description: Represent a wavefront as a sum of Legendre polynomials using
%              a matrix inversion.
% 
% This function attempts to solve the a_i's in equation,
% 
%                     M
%                     __
%                    \
%  phi(rho,theta) =  /__  a_i * P_i(x,y)
%                    i=1
% 
% where the P_i(x,y)'s are the Legendre polynomials from function 
% file, phi is the wavefront to be represented as a sum of Legendre 
% polynomials, the a_i's are the legendre coefficients, and M is the number
% of legendrepolynomials to use. (up to 7)
%
% Input:    phi - Phase to be represented as a sum of  legendrepolynomials
%                 that must be an nXmn array (rectangular)
%           (optional) M - Number of legendre polynomials to use (Default = 12)
% Output:   a - legendre coefficients (a_i's) as a vector
% S.Berujon 
% detourne des coefs de zernike sscript

% use of orthonormal polynoms
if nargin == 1
    M = 7;
end
if M > 7
    error('20? Why so high?.');
end

    x = -1:1/(128-1/2):1;

    leg = cell(M,2);
      
    

        leg{1,1} = repmat(ones(1,length(x)),256,1);   
        leg{1,2} = leg{1,1}';
        
        leg{2,1} = repmat((3^(1/2)).*x,256,1);   
        leg{2,2} = leg{2,1}';
        
        leg{3,1} = repmat(((5^(1/2))/2).*(3*(x.^2)-1),256,1);   
        leg{3,2} = leg{3,1}';
        
        leg{4,1} = repmat(((7^(1/2))/2).*(5*(x.^3)-3.*x),256,1);   
        leg{4,2} = leg{4,1}';
        
        leg{5,1} = repmat((3/8).*(35*(x.^4)-30.*(x.^2)+3),256,1);   
        leg{5,2} = leg{5,1}';
        
        leg{6,1} = repmat((11^(1/2))/8.*(63*(x.^5)-70.*(x.^3)+15.*x),256,1);   
        leg{6,2} = leg{6,1}';
        
        leg{7,1} = repmat((13^(1/2))/16.*(231*(x.^6)-315.*(x.^4)+105.*(x.^2)-5),256,1);   
        leg{7,2} = leg{7,1}';
        
        
        phi_size = size(phi);
    
        phi = reshape(phi,phi_size(1)*phi_size(2),1);
        Z = nan(phi_size(1)*phi_size(2),7);
        
        
        for k=1:7
            Z(:,k) = reshape(imresize(leg{k,1},[phi_size(1) phi_size(2)]),phi_size(1)*phi_size(2),1);
        end
        a_v = pinv(Z)*phi;
        
        
        for k=1:7
            Z(:,k) = reshape(imresize(leg{k,2},[phi_size(1) phi_size(2)]),phi_size(1)*phi_size(2),1);
        end
        a_h = pinv(Z)*phi;
        
    
        a=[a_v a_h];
        
reconst = zeros(phi_size(1),phi_size(2));
        
        
for kk = 1: M 
    reconst = a_v(kk) *imresize(leg{kk,1},[phi_size(1) phi_size(2)]) + reconst; 
end;

        
for kk = 1: M 
    reconst = a_h(kk) *imresize(leg{kk,2},[phi_size(1) phi_size(2)]) + reconst; 
end;        
        
        