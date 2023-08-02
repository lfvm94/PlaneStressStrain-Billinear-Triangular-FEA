function [es, st, stressXX_nodes_Elem, stressYY_nodes_Elem,stressXY_nodes_Elem,...
        strainXX_nodes_Elem,strainYY_nodes_Elem,strainXY_nodes_Elem]=...
        StressStrainCompatibilityTriang3(nElements,Eds,ngp,E,v,typeanalysis,...
        ex,ey,Edof,Coord)
%-------------------------------------------------------------------
% [stress,strain,stress_nodes,strain_nodes ] = ...
%  triang3StressStrain( ex, ey, ep, ue, ngp)
%-------------------------------------------------------------------
% PURPOSE
%  Compute the strain and stress at the center AND at each node
%  of a three-node-triangular finite element with gauss integration
%
% INPUT:  ex = [ x1 x2 x3 ]     element nodal x-coordinates
%         ey = [ y1 y2 y3 ]     element nodal y-coordinates
%
%         Eds:                  nodal displacements per node
%       
%         E:                    Modulus of elasticity
%         v:                    Poisson's ratio
% 
%         typeanalysis          Type of analysis: 1 - plane stress
%                                                 2 - plane strain
%
%         ngp                   Number of Gauss Points. Only 1 is
%                               allowed for this 2D linear triangle
%
%         Edof:                 Topology matrix: size: nelem x 7
%
%         Coord:                node coordinates [x,y]
%
%------------------------------------------------------------------
% OUTPUT: es : stress tensor at the center of element
%         st : strain tensor at the center of element
%
%         stressXX_nodes_Elem,  
%         stressYY_nodes_Elem,
%         stressXY_nodes_Elem:  stresses at each node of element
%
%         strainXX_nodes_Elem,  
%         strainYY_nodes_Elem,
%         strainXY_nodes_Elem:  strains at each node of element
%------------------------------------------------------------------
%
% CREATED by Luis F. Verduzco (2023-07-31)
%------------------------------------------------------------------

es=zeros(nElements,3);
st=zeros(nElements,3);

stressXX_nodes_Elem=zeros(nElements,3);
stressYY_nodes_Elem=zeros(nElements,3);
stressXY_nodes_Elem=zeros(nElements,3);

strainXX_nodes_Elem=zeros(nElements,3);
strainYY_nodes_Elem=zeros(nElements,3);
strainXY_nodes_Elem=zeros(nElements,3);

ngaussp=ngp;
Node_Stress=zeros(length(Coord(:,1)),3);
Node_Strain=zeros(length(Coord(:,1)),3);
Node_NoElem=zeros(length(Coord(:,1)),1);
for i=1:nElements
    ue=Eds(i,:);
    eps = [ typeanalysis E v ];

    [stress, strain, stress_nodes, strain_nodes] = triang3StressStrain...
        (ex(i,:),ey(i,:),eps,ue,ngaussp);

    stressXX_nodes_Elem(i,:)=stress_nodes(1,:);
    stressYY_nodes_Elem(i,:)=stress_nodes(2,:);
    stressXY_nodes_Elem(i,:)=stress_nodes(3,:);
    
    strainXX_nodes_Elem(i,:)=strain_nodes(1,:);
    strainYY_nodes_Elem(i,:)=strain_nodes(2,:);
    strainXY_nodes_Elem(i,:)=strain_nodes(3,:);
    
    % Sum stresses on each node of each element
    for j=1:3
        % Stresses XX
        Node_Stress(Edof(i,1 + 2*j)/2,1)=Node_Stress(Edof(i,1+2*j)/2,1)+...
                                       stressXX_nodes_Elem(i,j);
        % Stresses YY
        Node_Stress(Edof(i,1+2*j)/2,2)=Node_Stress(Edof(i,1+2*j)/2,2)+...
                                       stressYY_nodes_Elem(i,j);
        % Stresses XY
        Node_Stress(Edof(i,1+2*j)/2,3)=Node_Stress(Edof(i,1+2*j)/2,3)+...
                                       stressXY_nodes_Elem(i,j);
        % Strains XX
        Node_Strain(Edof(i,1 + 2*j)/2,1)=Node_Strain(Edof(i,1+2*j)/2,1)+...
                                       strainXX_nodes_Elem(i,j);
        % Strains YY
        Node_Strain(Edof(i,1+2*j)/2,2)=Node_Strain(Edof(i,1+2*j)/2,2)+...
                                       strainYY_nodes_Elem(i,j);
        % Strains XY
        Node_Strain(Edof(i,1+2*j)/2,3)=Node_Strain(Edof(i,1+2*j)/2,3)+...
                                       strainXY_nodes_Elem(i,j);
                                   
        % To count the times that a stress value is added to each node
        % Note: this data will be used to compute the average stress on
        % each node of the mesh
        Node_NoElem(Edof(i,1+2*j)/2)=Node_NoElem(Edof(i,1+2*j)/2)+1;
    end
    
    % Stresses and strains at the center of each element
    es(i,:)=stress;
    st(i,:)=strain;
end

% Computing the average stress on each node for compatibility
Node_Stress=Node_Stress./Node_NoElem(:,1);
Node_Strain=Node_Strain./Node_NoElem(:,1);

% Returning the new average nodal values to the stress arrays
for i=1:nElements
    for j=1:3
        node=Edof(i,1+2*j)/2;
        stressXX_nodes_Elem(i,j)=Node_Stress(node,1);
        stressYY_nodes_Elem(i,j)=Node_Stress(node,2);
        stressXY_nodes_Elem(i,j)=Node_Stress(node,3);
        
        strainXX_nodes_Elem(i,j)=Node_Strain(node,1);
        strainYY_nodes_Elem(i,j)=Node_Strain(node,2);
        strainXY_nodes_Elem(i,j)=Node_Strain(node,3);
    end
end
