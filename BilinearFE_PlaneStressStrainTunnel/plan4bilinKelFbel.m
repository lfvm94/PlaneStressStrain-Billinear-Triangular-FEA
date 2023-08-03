function [ Kel, fe ] = plan4bilinKelFbel( ex, ey, ep, eq )
% [ Ke, fe ] = plan4bilinKelFbel( ex, ey, ep, eq )
%-------------------------------------------------------------
% PURPOSE
%  Compute the stiffness matrix and element external force vector
%  for a bilinear plane stress or plane strain element including 
%  influence of out-of-plane stress (or strain)
%
% INPUT:  ex = [ x1 x2 x3 x4 ]         element nodal x-coordinates
%         ey = [ y1 y2 y3 y4 ]         element nodal y-coordinates
%
%         ep = [ ptype t ngp Emod ny ] ptype = 1 for plane stress = 1, or
%                                      ptype = 2 for plane strain
%                                        t: thickness
%                                      ngp: number of Gauss points in each
%                                           direction ( ksi and eta)
%                                           (ngp = 1 or 2 or 3)
%                                     Emod: Modulus of elasticity
%                                       ny: Poisson's ratio
%
%         eq = [ bx;                bx: body force x-dir
%                by ]               by: body force y-dir
%--------------------------------------------------------------
% OUTPUT: Kel : element stiffness matrix (8 x 8)
%         fel : equivalent nodal forces (8 x 1)
%--------------------------------------------------------------
% 
%
% MODIFIED for MEF by Luis Verduzco 2021-09-11
%--------------------------------------------------------------
%
% Pick out parameters from ep to facilitate readability of the code

ptype = ep(1);
t     = ep(2);
ngp   = ep(3); % Total gauss points
E  = ep(4);
v    = ep(5);

if ptype==1 % Plane stress
    D=(E/(1-v^2))*[1 v 0;
                    v 1 0;
                    0 0 (1-v)/2];
elseif ptype==2 % Plane strain
    d=(1-v)/(1-2*v);
    b=v/(1-2*v);
    
    D=(E/(1+v))*[d b 0;
                 b d 0;
                 0 0 0.5];
end


%  Tolerance for the Jacobian determinant
minDetJ = 1.e-16;

ngp   = ngp^2; % Total gauss points = ( NoGaussPoits per direction )^2

if ngp == 4
 
    intWeight=[1,1;
                1,1];


    GaussPoints=[-0.57735026918962,-0.57735026918962;
                  0.57735026918962,0.57735026918962];
             
else
 
    error('Only the option of 2 Gauss Points in each direction is allowed')

end

pvc=eq(2);

fe=zeros(8,1);
Kel=zeros(8);
for punto_gaus_xsi=1:ngp^0.5
    for punto_gaus_eta=1:ngp^0.5
        
        % Compute derivatives (with respect to xsi and eta) of the
        % shape functions at coordinate (xsi,eta). Since the element is
        % isoparametic, these are also the derivatives of the basis functions.
        xsi=GaussPoints(punto_gaus_xsi,1);
        weightXsi=intWeight(punto_gaus_xsi,1);
        
        eta=GaussPoints(punto_gaus_eta,2);
        weightEta=intWeight(punto_gaus_eta,2); 
    
        Ne=[(xsi-1).*(eta-1)/4  -(1+xsi).*(eta-1)/4 ... 
            (xsi+1).*(eta+1)/4 -(xsi-1).*(1+eta)/4]; % 1x4

        dNr=[-(1-eta)/4 (1-eta)/4 (1+eta)/4 -(1+eta)/4; % d(xsi)/d(eta)
             -(1-xsi)/4 -(1+xsi)/4 (1+xsi)/4 (1-xsi)/4]; % d(eta)/d(xsi)
     
        %  Use shape function derivatives and element vertex coordinates 
        %  to establish the Jacobian matrix.
    
        JT=dNr*[ex;ey]';

        detJ=det(JT);
        if detJ < minDetJ
            disp('Error. The Jacobian matrix is close to singular')
            break;
        end
        % Determinant seems OK - invert the transpose of the Jacobian
        JTinv=inv(JT);
        
        % Compute derivatives with respect to x and y, of all basis functions,
        dNxy=JTinv*dNr;

        Bxsieta=[dNxy(1,1) 0 dNxy(1,2) 0 dNxy(1,3) 0 dNxy(1,4) 0;
                0 dNxy(2,1) 0 dNxy(2,2) 0 dNxy(2,3) 0 dNxy(2,4);
                dNxy(2,1) dNxy(1,1) dNxy(2,2) dNxy(1,2) dNxy(2,3) dNxy(1,3) dNxy(2,4) dNxy(1,4)];

        % Compute the contribution to element stiffness and load matrix 
        % from current Gauss point
        Kel=Kel+Bxsieta'*D*Bxsieta*detJ*weightXsi*weightEta*t;

        % Gravity forces
        fe(2:2:8)=fe(2:2:8)+Ne'*detJ*weightXsi*weightEta*t*pvc;
    end
end
