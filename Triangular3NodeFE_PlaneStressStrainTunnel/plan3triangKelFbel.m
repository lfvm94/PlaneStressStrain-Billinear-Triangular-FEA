function [ Kel, fe ] = plan3triangKelFbel( ex, ey, ep, eq )
%--------------------------------------------------------------
% [ Ke, fe ] = plan3triangKelFbel( ex, ey, ep, eq )
%--------------------------------------------------------------
% PURPOSE
%  Compute the stiffness matrix and element external force vector
%  for a three node triangular plane stress or plane strain element.
%
% INPUT:  ex = [ x1 x2 x3 ]         element nodal x-coordinates
%         ey = [ y1 y2 y3]         element nodal y-coordinates
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
% OUTPUT: Kel : element stiffness matrix (6 x 6)
%         fel : equivalent nodal forces (6 x 1)
%--------------------------------------------------------------
%
% MODIFIED for MEF by Luis Verduzco 2021-09-11
%--------------------------------------------------------------
%
% Pick out parameters from ep to facilitate readability of the code

ptype = ep(1);
t     = ep(2);
ngp   = ep(3); % Total gauss points. For this case, only 1 is allowed
E  = ep(4);
v    = ep(5);

%% Compute the constitutive matrix D
if ptype==1
    D=(E/(1-v^2))*[1 v 0;
                    v 1 0;
                    0 0 (1-v)/2];
elseif ptype==2
    d=(1-v)/(1-2*v);
    b=v/(1-2*v);
    
    D=(E/(1+v))*[d b 0;
                 b d 0;
                 0 0 0.5];
end


%  Tolerance for the Jacobian determinant
minDetJ = 1.e-16;

if ngp == 1
 
    intWeight=[0.5,0.5];
    GaussPoints=[1/3,1/3];
   
else
 
    error('Only 1 Gauss Point applies')

end

pvc=eq(2);

%  Initialize Ke and fe with zeros for all of their elements
fe=zeros(6,1);
Kel=zeros(6);
for gaussPointXsiEta=1:ngp

    % Compute derivatives (with respect to xsi and eta) of the shape
    % functions at coordinate (xsi,eta). Since the element is isoparametic,
    % these are also the derivatives of the basic functions.
    xsi=GaussPoints(gaussPointXsiEta,1);
    weightXsi=intWeight(gaussPointXsiEta,1);

    eta=GaussPoints(gaussPointXsiEta,2);
    weightEta=intWeight(gaussPointXsiEta,2); 

    Ne=[1-xsi-eta xsi eta]; % 1x3

    dNr=[-1 1 0; % d(xsi)/d(eta)
         -1 0 1];  % d(eta)/d(xsi)

    %  Use shape function derivatives and element vertex coordinates 
    %  to compute the Jacobian matrix.

    JT=dNr*[ex;ey]';

    detJ=det(JT);
    if detJ < minDetJ
        disp('Error. The Jacobian matrix is close to singular !');
        return;
    end
    
    % Determinant seems OK - invert the transpose of the Jacobian
    JTinv=inv(JT);
    
    % Compute derivatives with respect to x and y, of all basis functions,
    dNxy=JTinv*dNr;

    Bxsieta=[dNxy(1,1) 0 dNxy(1,2) 0 dNxy(1,3) 0;
            0 dNxy(2,1) 0 dNxy(2,2) 0 dNxy(2,3);
            dNxy(2,1) dNxy(1,1) dNxy(2,2) dNxy(1,2) dNxy(2,3) dNxy(1,3)];

    % Compute the contribution to element stiffness and load matrix 
    % from current Gauss point
    Kel=Kel+Bxsieta'*D*Bxsieta*detJ*weightXsi*weightEta*t;

    % Gravity forces
    fe(2:2:6)=fe(2:2:6)+Ne'*detJ*weightXsi*weightEta*t*pvc;
    
end
