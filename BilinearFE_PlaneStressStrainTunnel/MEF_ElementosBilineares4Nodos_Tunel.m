clear all
clc

%% Tunnel mesh and geometry
H=1000; % Tunnel's height (cm)
B=800; % Half tunnel's width
Z=5000; % Depth at which the upper tunnel's boundary is under water (cm)
z=0.1;
D=300; % Thickness of tunnel's pavement
b=400; % Interior tunnel's width
h=400; % Interior tunnel's height
rs=40; % Interior tunnel's edge radius

[Edof,Coord,Ex,Ey,LeftSide_nodes,TopSide_nodes,RightSide_nodes,...
    BottomSide_nodes]=TunnelMeshGen(H,B,D,b,h,rs,30,50,2);

t=1; % constant thickness (cm)
nElements=length(Edof(:,1)); % Number of Finite Elements
nnodes=length(Coord); % Number of nodes of the mesh
ndof=2*nnodes;

%% Materials
fc=300; % concrete's compressive strength (Kg/cm2) 
E=11000*sqrt(fc); % Modulus of Elasticity of Concrete
v=0.2; % Poisson ratio
pvc=-2400e-6; % volumetric weight 

%% Type of analysis: plane stress (1) or plane strain (2)

typeAnalysis=1; % 1: plane stress, 2: plane strain

eq=[ 0; pvc]; % self weight loads on each global axis direction:
            % [x(horizontal), y(gravity)]
            
%% Gauss points for analysis
ngp=2; % Number of Gauss points for numerical integration for each
       % Billinear Finite Element. Only 2 gauss points are allowed 
       % for this element so that the integration is best approximated

ep=[typeAnalysis t ngp E v ];

K=zeros(ndof);
f=zeros(ndof,1);

for i=1:nElements
        
    [ Ke, fe ] = plan4bilinKelFbel(Ex(i,:),Ey(i,:),ep,eq);
    
     for j=2:9
         
         K(Edof(i,j),Edof(i,2))=K(Edof(i,j),Edof(i,2))+Ke(j-1,1);
         K(Edof(i,j),Edof(i,3))=K(Edof(i,j),Edof(i,3))+Ke(j-1,2);
         K(Edof(i,j),Edof(i,4))=K(Edof(i,j),Edof(i,4))+Ke(j-1,3);
         K(Edof(i,j),Edof(i,5))=K(Edof(i,j),Edof(i,5))+Ke(j-1,4);
         K(Edof(i,j),Edof(i,6))=K(Edof(i,j),Edof(i,6))+Ke(j-1,5);
         K(Edof(i,j),Edof(i,7))=K(Edof(i,j),Edof(i,7))+Ke(j-1,6);
         K(Edof(i,j),Edof(i,8))=K(Edof(i,j),Edof(i,8))+Ke(j-1,7);
         K(Edof(i,j),Edof(i,9))=K(Edof(i,j),Edof(i,9))+Ke(j-1,8);

     end

     f(Edof(i,2))=f(Edof(i,2))+fe(1);
     f(Edof(i,3))=f(Edof(i,3))+fe(2);
     f(Edof(i,4))=f(Edof(i,4))+fe(3);
     f(Edof(i,5))=f(Edof(i,5))+fe(4);
     f(Edof(i,6))=f(Edof(i,6))+fe(5);
     f(Edof(i,7))=f(Edof(i,7))+fe(6);
     f(Edof(i,8))=f(Edof(i,8))+fe(7);
     f(Edof(i,9))=f(Edof(i,9))+fe(8);
end


%% SOLVER
%% Restrictions (constraints)
%   [        DOF,           Prescribed displacement      ]
bc=[BottomSide_nodes*2-1 zeros(length(BottomSide_nodes),1);
    BottomSide_nodes*2   zeros(length(BottomSide_nodes),1);
    LeftSide_nodes*2-1   zeros(length(LeftSide_nodes),1)];

%% External forces
% Distribution of forces at the mesh's boundaries:
fwater=zeros(ndof,1);

% Upper boundary
ww=1e-3; % Unit weight of water (Kg/cm3)
wub=Z*ww; % kg/cm2 (uniformly distributed load magnitude at the top)

for i=1:length(TopSide_nodes)-1
    % Computing the distance between each node at the boundary in question
    x1=Coord(TopSide_nodes(i),1);
    x2=Coord(TopSide_nodes(i+1),1);
    
    y1=Coord(TopSide_nodes(i),2);
    y2=Coord(TopSide_nodes(i+1),2);
    
    lenBoundTop=((x1-x2)^2+(y1-y2)^2)^0.5;
    
    % Apply distributed loads at the boundary's nodes
    fed=[0;
        -wub*t*lenBoundTop/2;
        0;
        -wub*t*lenBoundTop/2]; % only gravity loads are applied
                                     % for this case (negative Y direction)

    fwater(2*TopSide_nodes(i))=fwater(2*TopSide_nodes(i))+...
                               fed(2);
                           
    fwater(2*TopSide_nodes(i+1))=fwater(2*TopSide_nodes(i+1))+...
                                 fed(4);
end

% Right boundary
for i=1:length(RightSide_nodes)-1
    % Computing the distance between each node at the boundary in question
    x1=Coord(RightSide_nodes(i),1);
    x2=Coord(RightSide_nodes(i+1),1);
    
    y1=Coord(RightSide_nodes(i),2);
    y2=Coord(RightSide_nodes(i+1),2);
    
    lenBoundRight=((x1-x2)^2+(y1-y2)^2)^0.5;

    % Apply distributed loads at the boundary's nodes
    za=(Z+(H-(y1+y2)*0.5)); % Average finite boundary depth
    wrb=za*ww;
    
    fed=[-wrb*t*lenBoundRight/2;
        0;
        -wrb*t*lenBoundRight/2; % forces at applied only for the
        0];                             % horizontal direction due to the
                                        % lateral water pressure

    fwater(2*RightSide_nodes(i)-1)=fwater(2*RightSide_nodes(i)-1)+...
                                   fed(1);
    fwater(2*RightSide_nodes(i+1)-1)=fwater(2*RightSide_nodes(i+1)-1)+...
                                     fed(3);

end

% Summ external forces and internal forces
f=f+fwater;

% Solving the system of equations
[d,r]=solveq(K,f,bc); % This is a CALFEM's function
                      % To download CALFEM go to its repository:
                      % https://github.com/CALFEM/calfem-matlab
                      
%% Graphic solutions 

%-----Undeformed mesh-----%
figure(1);
xlabel('Width [cm]')
ylabel('Height [cm]')
plotpar=[1 1 0];
eldraw2(Ex,Ey,plotpar); % This is a CALFEM's function
                      % To download CALFEM go to its repository:
                      % https://github.com/CALFEM/calfem-matlab

%----Deformed mesh----%
figure(1)
sfac=1000;
Ed=extract(Edof,d);
plotpar=[1 3 2];
eldisp2(Ex,Ey,Ed,plotpar,sfac); % This is a CALFEM's function
                      % To download CALFEM go to its repository:
                      % https://github.com/CALFEM/calfem-matlab
title(strcat('Deformed - Undeformed mesh', 'Scale x ',num2str(sfac)));

%% Stresses [Kg/cm2] & Strains 
[es, st, stressXX_nodes_Elem, stressYY_nodes_Elem,stressXY_nodes_Elem,...
strainXX_nodes_Elem,strainYY_nodes_Elem,strainXY_nodes_Elem]=...
StressStrainBillinear4(nElements,Ed,ngp,E,v,typeAnalysis,Ex,Ey,Edof,...
Coord);

Edu=extract(Edof,d) ;

%%%-----------------------------------------------------------------------
%--------------------Displacements distribution plot--------------------%
%%%-----------------------------------------------------------------------

figure(2);
fill(Ex',Ey',Edu(:,[1 3 5 7])');
axis equal;  
axis image;
shading interp;
colormap 'jet';
set(gca,'fontsize',13);
xlabel('Width [cm]')
ylabel('Height [cm]')
title('Displacement xx Distribution');

figure(3);
fill(Ex',Ey',Edu(:,[2 4 6 8])');
axis equal;  
axis image;
shading interp;
colormap 'jet';
set(gca,'fontsize',13);
xlabel('Width [cm]')
ylabel('Height [cm]')
title('Displacement yy Distribution');

%%               Stresses at the center of each element
%%%-----------------------------------------------------------------------
%---------------------------Stresses Plots-------------------------------%
%%%-----------------------------------------------------------------------

%%%--------------------Plane stresses xx direction---------------------%%%
figure(4);
fill(Ex',Ey',es(:,1));
axis equal;
axis image;
shading flat;
colormap 'jet';
xlabel('Width [cm]')
ylabel('Height [cm]')
title(strcat('Plane stresses xx distribution (Kg/cm2) - ',...
             'Stress at the center of elements'))

%%%-----------------------------------------------------------------------
%%%-----------------------------------------------------------------------
%%%--------------------Plane stresses yy direction---------------------%%%

figure(5);
fill(Ex',Ey',es(:,2));
axis equal;
axis image;
shading flat;
colormap 'jet';
xlabel('Width [cm]')
ylabel('Height [cm]')
title(strcat('Plane stresses yy distribution (Kg/cm2) - ',...
             'Stress at the center of elements'))
hold on

%%%--------------------Plane stresses xy direction---------------------%%%

figure(6);
fill(Ex',Ey',es(:,3));
axis equal;
axis image;
shading flat;
colormap 'jet';
xlabel('Width [cm]')
ylabel('Height [cm]')
title(strcat('Plane stresses xy distribution (Kg/cm2) - ',...
             'Stress at the center of elements'))
hold on

%%             Stresses interpolated on each element's nodes
%%%------------------Plots with interpolated contours--------------------
%%%-----------------------------------------------------------------------

figure(9);
fill(Ex',Ey',stressXX_nodes_Elem');
axis equal;
axis image;
shading interp;
colormap 'jet';
set(gca,'fontsize',13);
xlabel('Width [cm]')
ylabel('Height [cm]')
title('Plane stresses xx distribution (Kg/cm2)')
hold on
   
figure(10);
fill(Ex',Ey',stressYY_nodes_Elem');
axis equal;
axis image;
shading interp;
colormap 'jet';
set(gca,'fontsize',13);
xlabel('Width [cm]')
ylabel('Height [cm]')
title('Plane stresses yy distribution (Kg/cm2)')
hold on
   
figure(11);
fill(Ex',Ey',stressXY_nodes_Elem');
axis equal;
axis image;
shading interp;
colormap 'jet';
set(gca,'fontsize',13);
xlabel('Width [cm]')
ylabel('Height [cm]')
title('Plane stresses xy distribution (Kg/cm2)')
hold on
   
