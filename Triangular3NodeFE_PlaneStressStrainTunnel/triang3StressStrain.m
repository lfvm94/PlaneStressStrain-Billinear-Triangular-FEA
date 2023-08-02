  function [stress,strain,stress_nodes,strain_nodes ] = ...
            triang3StressStrain( ex, ey, ep, ue, ngp)
%-------------------------------------------------------------------
% [stress,strain,stress_nodes,strain_nodes ] = ...
%  triang3StressStrain( ex, ey, ep, ue, ngp)
%-------------------------------------------------------------------
% PURPOSE
%  Compute the strain and stress at the center and at each node
%  of a three node triangular finite element with gauss integration
%
% INPUT:  ex = [ x1 x2 x3 ]    element nodal x-coordinates
%         ey = [ y1 y2 y3 ]    element nodal y-coordinates
%
%         ep = [ ptype Emod ny ]  ptype = 1 for plane stress = 1, or
%                                 ptype = 2 for plane strain                                 
%                                 Emod: Modulus of elasticity
%                                   ny: Poisson's ratio
%
%         u_e                     element displacement matrix,
%                                 size = ( 1, 6 )
%
%         ngp                     Number of Gauss Points. Only 1 is
%                                 allowed for this 2D linear triangle
%------------------------------------------------------------------
% OUTPUT: stress : stress tensor at the center of element
%         strain : strain tensor at the center of element
%         stress_nodes: stresses at each node of element
%         strain_nodes: strains at each node of element
%------------------------------------------------------------------
% 
% MODIFIED by Dimosthenis Floros (2015-12-01)
% MODIFIED by Luis F. Verduzco (2023-07-31)
%------------------------------------------------------------------
%
% Pick out parameters from ep to facilitate readability of the code

ptype=ep(1);
E=ep(2);
v= ep(3);

% Determine constitutive matrix D for plane strain or plane stress

if ptype == 1
    D=E/(1-v^2)* [1  v   0;
                  v  1   0;
                  0  0 (1-v)/2];
 
else
    D=E/(1+v)/(1-2*v)* [1-v  v    v;
                         v  1-v   v;
                         0    0   (1-2*v)/2];

end

% For each number of gauss point the coordinates (xsi, eta) in the 
% parent domain are:

if ngp == 1
 
    intWeight=[0.5,0.5];
    GaussPoints=[1/3,1/3];
               
else
 
    error('Only 1 Gauss Point applies')

end
     
st=zeros(3,ngp);
str=zeros(3,ngp);
stress_nodes=zeros(3,3); % to store the nodal stresses of each element
strain_nodes=zeros(3,3); % to store the nodal strains of each element
for gauss_point=1:ngp

    xsi =GaussPoints(gauss_point,1);
    eta =GaussPoints(gauss_point,2);

    % Compute the derivatives (with respect to xsi and eta) of the
    % shape functions at coordinate (xsi,eta)

    dNr=[-1 1 0; % d(xsi)/d(eta)
         -1 0 1]; % d(eta)/d(xsi)
    
    % Compute Jacobian matrix and invert the transpose of the Jacobian

    JT=dNr*[ex;ey]';

    detJ=det(JT);
    JTinv=inv(JT);

    % Compute derivatives with respect to x and y, of all basis functions

    dNxy=JTinv*dNr;

    % Use the derivatives of the shape functions to compute the element
    % B-matrix, Be

    Bxsieta=[dNxy(1,1) 0 dNxy(1,2) 0 dNxy(1,3) 0;
             0 dNxy(2,1) 0 dNxy(2,2) 0 dNxy(2,3);
             dNxy(2,1) dNxy(1,1) dNxy(2,2) dNxy(1,2) dNxy(2,3) dNxy(1,3)];

    % Compute the strain and store it in strain

    str(:,gauss_point) = Bxsieta*ue';
    st(:,gauss_point) = D*str(:,gauss_point);
    
    % One loop for each node
    for i=1:3
        xsi_n=1/xsi*ex(i);
        eta_n=1/eta*ey(i);

        Ne=[1-xsi_n-eta_n xsi_n eta_n]; % 1 x 3 

        % Compute the stress and strains
        esf=Ne.*st(:,1);
        def=Ne.*str(:,1);

        stress_nodes(:,i)=sum(esf')';
        strain_nodes(:,i)=sum(def');
    end
end

% Computation of stresses and strains at the center of each element
% as the average of the values computed fo reach gauss point. This values
% are used when the parameter 'flat' is used on the plots

strains=zeros(3,1);
strains(1)=sum(str(1,:))/ngp;
strains(2)=sum(str(2,:))/ngp;
strains(3)=sum(str(3,:))/ngp;

strain=strains';

stresses=zeros(3,1);
stresses(1)=sum(st(1,:))/ngp;
stresses(2)=sum(st(2,:))/ngp;
stresses(3)=sum(st(3,:))/ngp;

stress=stresses';

% -------------------------------- End --------------------------------- %
