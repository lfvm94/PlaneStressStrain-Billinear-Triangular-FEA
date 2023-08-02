  function [stress,strain,stress_nodes,strain_nodes ] = ...
            plan4bilinStressStrain( ex, ey, ep, ue, ngp)
        
% [stress,strain,stress_nodes,strain_nodes ] = ...
%            plan4bilinStressStrain( ex, ey, ep, ue, ngp)
%-------------------------------------------------------------
% PURPOSE
%  Compute the strain and stress in the center of a quadrilateral
%  bilear finite element
%
% INPUT:  ex = [ x1 x2 x3 x4 ]    element nodal x-coordinates
%         ey = [ y1 y2 y3 y4 ]    element nodal y-coordinates
%
%         ep = [ ptype Emod ny ]  ptype = 1 for plane stress = 1, or
%                                 ptype = 2 for plane strain                                 
%                                 Emod: Modulus of elasticity
%                                   ny: Poisson's ratio
%
%         ue                      element displacement matrix,
%                                 size = ( 1, 8 )
%
%         ngp                     Number of Gauss Points from 1 or 2
%-------------------------------------------------------------
% OUTPUT: stress : stress tensor at the center of element
%         strain : strain tensor at the center of element
%         stress_nodes: stresses at each node of element
%         strain_nodes: strains at each node of element
%--------------------------------------------------------------
% 
%
% MODIFIED by Dimosthenis Floros (20151201)
% MODIFIED by Luis F. Verduzco (20210910)
%--------------------------------------------------------------
%
% Pick out parameters from ep to facilitate readability of the code

ptype=ep(1);
E=ep(2);
v= ep(3);

% Determine constitutive matrix D for plane strain or plane stress

if ptype == 1
    D=E/(1-v^2)*[1  v   0;
                v   1   0;
                0   0 (1-v)/2];
 
else
    D=E/(1+v)/(1-2*v)*[1-v  v   v;
                       v  1-v   v;
                       0    0   (1-2*v)/2];

end


% For 1 gauss point the coordinates (xsi, eta) in the parent domain are:


ngp=ngp^2;

if ngp == 4
             % xsi eta
    intWeight=[1,1;
               1,1];

                 % xsi                 eta
    GaussPoints=[-0.57735026918962, -0.57735026918962;
                 0.57735026918962,  0.57735026918962];
         
else
 
    error('Only the option of 2 Gauss Points in each direction apply')

end

gauss_point=0; % counter
               % 2 4
               % 1 3
st=zeros(3,ngp);
str=zeros(3,ngp);
stress_nodes=zeros(3,4); % esfuerzos en nodos de elemento
strain_nodes=zeros(3,4); % strains en nodos de elemento
for gauspoint_xsi=1:ngp^0.5
    
    for gauspoint_eta=1:ngp^0.5
        gauss_point=gauss_point+1;
        
        xsi =GaussPoints(gauspoint_xsi,1)+1e-9;
        eta =GaussPoints(gauspoint_eta,2)+1e-9;
    
        % Compute the derivatives (with respect to xsi and eta) of the
        % shape functions at coordinate (xsi,eta).

        dNr=[-(1-eta)/4 (1-eta)/4 (1+eta)/4 -(1+eta)/4;
             -(1-xsi)/4 -(1+xsi)/4 (1+xsi)/4 (1-xsi)/4];

        % Compute Jacobian matrix and invert the transpose of the Jacobian

        JT=dNr*[ex;ey]';

        detJ=det(JT);
        JTinv=inv(JT);

        % Compute derivatives with respect to x and y, of all basis functions

        dNxy=JTinv*dNr;

        % Use the derivatives of the shape functions to compute the element
        % B-matrix, Be

        Bxsieta=[dNxy(1,1) 0 dNxy(1,2) 0 dNxy(1,3) 0 dNxy(1,4) 0;
                 0 dNxy(2,1) 0 dNxy(2,2) 0 dNxy(2,3) 0 dNxy(2,4);
                 dNxy(2,1) dNxy(1,1) dNxy(2,2) dNxy(1,2) dNxy(2,3) dNxy(1,3) dNxy(2,4) dNxy(1,4)];

        % Compute the strain and store it in strain

        str(:,gauss_point) = Bxsieta*ue';
        st(:,gauss_point) = D*str(:,gauss_point);

        % Computing the stresses on each node:
        % For the case in which 2 gauss points are evaluated per axis:
        % 4 3
        % 1 2
        
        if gauss_point==1
            n=1;
            xsi_n=ex(n)/xsi;
            eta_n=ey(n)/eta;
            
            Ne=[(xsi_n-1).*(eta_n-1)/4  -(1+xsi_n).*(eta_n-1)/4 ... 
                (xsi_n+1).*(eta_n+1)/4 -(xsi_n-1).*(1+eta_n)/4]; % 1x4
         
            % Compute the stress and store it in stress
            esf=Ne.*st(:,gauss_point);
            def=Ne.*str(:,gauss_point);

            stress_nodes(:,n)=sum(esf')';
            strain_nodes(:,n)=sum(def');
            
        elseif gauss_point==2
            n=4;
            xsi_n=ex(n)/xsi;
            eta_n=ey(n)/eta;
            
            Ne=[(xsi_n-1).*(eta_n-1)/4  -(1+xsi_n).*(eta_n-1)/4 ... 
                (xsi_n+1).*(eta_n+1)/4 -(xsi_n-1).*(1+eta_n)/4]; % 1x4
         
            % Compute the stress and store it in stress
            esf=Ne.*st(:,gauss_point);
            def=Ne.*str(:,gauss_point);

            stress_nodes(:,n)=sum(esf')';
            strain_nodes(:,n)=sum(def');
            
        elseif gauss_point==3
            n=2;
            xsi_n=ex(n)/xsi;
            eta_n=ey(n)/eta;
            
            Ne=[(xsi_n-1).*(eta_n-1)/4  -(1+xsi_n).*(eta_n-1)/4 ... 
                (xsi_n+1).*(eta_n+1)/4 -(xsi_n-1).*(1+eta_n)/4]; % 1x4
         
            % Compute the stress and store it in stress
            esf=Ne.*st(:,gauss_point);
            def=Ne.*str(:,gauss_point);

            stress_nodes(:,n)=sum(esf')';
            strain_nodes(:,n)=sum(def');
            
        elseif gauss_point==4
            n=3;
            xsi_n=ex(n)/xsi;
            eta_n=ey(n)/eta;
            
            Ne=[(xsi_n-1).*(eta_n-1)/4  -(1+xsi_n).*(eta_n-1)/4 ... 
                (xsi_n+1).*(eta_n+1)/4 -(xsi_n-1).*(1+eta_n)/4]; % 1x4
         
            % Compute the stress and store it in stress
            esf=Ne.*st(:,gauss_point);
            def=Ne.*str(:,gauss_point);

            stress_nodes(:,n)=sum(esf')';
            strain_nodes(:,n)=sum(def');
        end
        
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
