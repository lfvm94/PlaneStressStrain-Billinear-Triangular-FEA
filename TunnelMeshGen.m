
function [Edof,Coord,Ex,Ey,LeftSide_nodes,TopSide_nodes,RightSide_nodes,...
    BottomSide_nodes]=TunnelMeshGen(H,B,D,b,h,r,Nr,Nt,elemtype)
% Generates a rectangular mesh of
% call:
% [Edof,Coord,Ex,Ey,LeftSide_nodes,TopSide_nodes,RightSide_nodes,...
% BottomSide_nodes,TopLeftElem_no,TopRighty_node,h]=RectangleMeshGen...
% (Lx,Ly,Nt,Nr,elemtype)
% -------------------------------------------------------------------------
% Input -------------------------------------------------------------------
%
% H is the domain size in y-direction
%
% B is the domain size in x-direction
%
% D is the distance between bottom of domain and tunnel floor
%
% b is the tunnel (half) width 
%
% h is the tunnel maximum height
%
% r is the radius of corner fillet
%
% Nr is the number of element side lengths from tunnel to domain top/bottom edge, 
% i.e. defines the number of elements along the top left and bottom left boundary
% (minimum 1)
%
% Nt is the number of element side lengths in tunnel inside
% (minimum 5)
%
% Set to 1 for CST elements(triangles) and set as 2 for 4-node bilinear
% elements
%
% Output ------------------------------------------------------------------
%
% Edof is the element topology matrix (for 2 dof per node)
%
% Coord is the global coordinate matrix
%
% Ex is the element nodal coordinates in x-direction
%
% Ey is the element nodal coordinates in y-direction
%
% LeftSide_nodes is a list of the nodes on the left boundary
% TopSide_nodes is a list of the nodes on the top boundary
% RightSide_nodes is a list of the nodes on the right boundary
% BottomSide_nodes is a list of the nodes on the bottom boundary
%
%--------------------------------------------------------------------------
% Developed for VSM167 FEM BASICS 2012 by Håkan Johansson
% Dept. Applied Mechanics, Chalmers University of Technology
%-------------------------------------------------------------

if H<0.9*(h+D) || B<0.9*b || b<r || h< r
    error('Inappropriate geometry chosen')
end

% divide tunnel boundary in suitable way ("divide largest element")
L_inside=b+(h-r)+pi*r/2 +(b-r);

%initial division
s_start=[0 b/L_inside (b+(h-r))/L_inside ...
        (b+(h-r)+pi*r/4)/L_inside (b+(h-r)+pi*r/2)/L_inside 1];
s=s_start;

% split largest edge-part in two
for i=6:Nt
   [~, ind]=max(diff(s));
   s=[s(1:ind) (s(ind)+s(ind+1))/2 s(ind+1:end)];
end

s_seg1=find(s<s_start(2));
s_seg2=find(s_start(2)<s & s< s_start(4));
s_seg3=find(s_start(4)< s);


X_squares=zeros(Nr+1,Nt+1);
Y_squares=zeros(Nr+1,Nt+1);

for i=s_seg1
   y_1=D;
   y_2=0;
   x_1=s(i)/s(length(s_seg1)+1)*b;
   x_2=s(i)/s(length(s_seg1)+1)*B;
   X_squares(:,i)=linspace(x_1,x_2,Nr+1);
   Y_squares(:,i)=linspace(y_1,y_2,Nr+1);
end
X_squares(:,s_seg1(end)+1)=linspace(b,B,Nr+1);
Y_squares(:,s_seg1(end)+1)=linspace(D,0,Nr+1);

for i=s_seg2
   x_2=B;
   y_2=(s(i)-s_start(2))/(s_start(4)-s_start(2))*H;
   if s(i)<s_start(3)
       x_1=b;
       y_1=(s(i)-s_start(2))/(s_start(3)-s_start(2))*(h-r)+D;
   else
      v=(s(i)-s_start(3))*L_inside/r;
      x_1=b-r+r*cos(v);
      y_1=(D+h)-r+r*sin(v);
   end
   
   X_squares(:,i)=linspace(x_1,x_2,Nr+1);
   Y_squares(:,i)=linspace(y_1,y_2,Nr+1);
end

X_squares(:,s_seg2(end)+1)=linspace(b-(1-1/sqrt(2))*r,B,Nr+1);
Y_squares(:,s_seg2(end)+1)=linspace(D+h-(1-1/sqrt(2))*r,H,Nr+1);

for i=s_seg3
   
   y_2=H;
   x_2=(1-(s(i)-s_start(4))/(1-s_start(4)))*B;
   if s(i)<s_start(5)
     v=(s(i)-s_start(3))*L_inside/r;
      x_1=b-r+r*cos(v);
      y_1=(D+h)-r+r*sin(v);
   else
       y_1=D+h;
       x_1=(1-(s(i)-s_start(5))/(1-s_start(5)))*(b-r);
   end
   
   X_squares(:,i)=linspace(x_1,x_2,Nr+1);
   Y_squares(:,i)=linspace(y_1,y_2,Nr+1);
end

%x_div=linspace(0,B,Nt+1);  % division in "squares" in x-direction

%y_div=linspace(0,H,Nt+1);

%[X_squares,Y_squares]=meshgrid(x_div,fliplr(y_div));

if elemtype==1      %triangles
    
    % introduce mid-nodes in each "square" element
    X_mids=(X_squares(1:end-1,1:end-1)+X_squares(1:end-1,2:end)+...
            X_squares(2:end,1:end-1)+X_squares(2:end,2:end))/4;
    Y_mids=(Y_squares(1:end-1,1:end-1)+Y_squares(1:end-1,2:end)+...
            Y_squares(2:end,1:end-1)+Y_squares(2:end,2:end))/4;
    
    X_coord=X_squares(1,:)';
    Y_coord=Y_squares(1,:)';
    for k=1:Nr   % each row of "squares"
        X_coord=[X_coord;X_mids(k,:)';X_squares(k+1,:)'];
        Y_coord=[Y_coord;Y_mids(k,:)';Y_squares(k+1,:)'];
    end
    
    
    Row_nod_top=[1:Nt;       % Each square split in 4 elements (layout "X")
        Nt+2:2*Nt+1;         % the top elements in row
        2:Nt+1]';
    
    Row_nod_left=[1:Nt;      % Left elements
        2*Nt+2:3*Nt+1;
        Nt+2:2*Nt+1]';
    
    Row_nod_right=[2:Nt+1;   % Right elements
        Nt+2:2*Nt+1;
        2*Nt+3:3*Nt+2]';
    
    Row_nod_bottom=[Nt+2:2*Nt+1;      % Bottom elements
        2*Nt+2:3*Nt+1;
        2*Nt+3:3*Nt+2]';
    
    Row_nod=[Row_nod_top;Row_nod_left;Row_nod_right;Row_nod_bottom]; % all elements in one row of "squares"
    
    Enod=[];
    parfor k=1:Nr   % each row of "squares"
        Enod=[Enod;Row_nod+(k-1)*(2*Nt+1)]; % next row has the same element/nodes, but are shifted 2*n+1 (=number of nodes on one row)
    end
    Enod=[(1:size(Enod,1))',Enod];   % add element numbers
    
    Edof=[Enod(:,1), 2*Enod(:,2)-1,2*Enod(:,2),2*Enod(:,3)-1,2*Enod(:,3),2*Enod(:,4)-1,2*Enod(:,4)];
    Ex=X_coord(Enod(:,2:end));
    Ey=Y_coord(Enod(:,2:end));
    
    % Characteristic element size is the mean of the longest distance between
    % two nodes in the element
    ElemDist=[sqrt((Ex(:,1)-Ex(:,2)).^2+(Ey(:,1)-Ey(:,2)).^2), ...
              sqrt((Ex(:,2)-Ex(:,3)).^2+(Ey(:,2)-Ey(:,3)).^2), ...
              sqrt((Ex(:,1)-Ex(:,3)).^2+(Ey(:,1)-Ey(:,3)).^2)];
  
    
elseif elemtype==2      % bilinear quadrilaterals
    
    X_coord=reshape(X_squares',[],1);
    Y_coord=reshape(Y_squares',[],1);
    
    Row_nod=[1:Nt;
        Nt+2:2*Nt+1;
        Nt+3:2*Nt+2;
        2:Nt+1]';
    
    Enod=[];
    parfor k=1:Nr   % each row of "squares"
        Enod=[Enod;Row_nod+(k-1)*(Nt+1)]; % next row has the same element/nodes, but are shifted 2*n+1 (=number of nodes on one row)
    end
    Enod=[(1:size(Enod,1))',Enod];   % add element numbers
    
    Edof=[Enod(:,1), 2*Enod(:,2)-1,2*Enod(:,2),2*Enod(:,3)-1,2*Enod(:,3),...
          2*Enod(:,4)-1,2*Enod(:,4),2*Enod(:,5)-1,2*Enod(:,5)];
    
    
    
    
    % Characteristic element size is the mean of the longest elements side
    % assume bilinear elements are rectangular
    Ex=X_coord(Enod(:,2:end));
    Ey=Y_coord(Enod(:,2:end));
    
    ElemDist=[sqrt((Ex(:,1)-Ex(:,2)).^2+(Ey(:,1)-Ey(:,2)).^2), ...
              sqrt((Ex(:,2)-Ex(:,3)).^2+(Ey(:,2)-Ey(:,3)).^2)];
   
    
    
end

Coord=[X_coord, Y_coord];

LeftSide_nodes=find(abs(X_coord-0)<100*eps);
RightSide_nodes=find(abs(X_coord-B)<100*eps);
BottomSide_nodes=find(abs(Y_coord-0)<100*eps);
TopSide_nodes=find(abs(Y_coord-H)<100*eps);


% To plot mesh
%eldraw2(Ex,Ey)

