function [es, st, stressXX_nodes_Elem, stressYY_nodes_Elem,...
    stressXY_nodes_Elem,strainXX_nodes_Elem,strainYY_nodes_Elem,...
    strainXY_nodes_Elem]=StressStrainBillinear4(nElements,Eds,ngp,...
    E,v,typeAnalysis,ex,ey,Edof,Coord)
%-------------------------------------------------------------------
% [es, st, stressXX_nodes_Elem, stressYY_nodes_Elem,...
%  stressXY_nodes_Elem,strainXX_nodes_Elem,strainYY_nodes_Elem,...
%  strainXY_nodes_Elem]=StressStrainBillinear4(nElements,Eds,ngp,...
%  E,v,typeAnalysis,ex,ey,Edof,Coord)
%-------------------------------------------------------------------
% PURPOSE
%  Compute the strain and stress at the at each node of a 
%  four-node-billinear finite element with gauss integration.
%
% INPUT:  ex = [ x1 x2 x3 x4]   element nodal x-coordinates
%         ey = [ y1 y2 y3 y4]   element nodal y-coordinates
%
%         Eds:                  nodal displacements per node
%       
%         E:                    Modulus of elasticity
%         v:                    Poisson's ratio
%
%         typeanalysis          Type of analysis: 1 - plane stress
%                                                 2 - plane strain
%
%         ngp                   Number of Gauss Points. Options are 1 
%                               or 2 for this 2D four-node billinear
%                               element
%
%         Edof:                 Topology matrix: size: nelem x 9
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

stressXX_nodes_Elem=zeros(nElements,4);
stressYY_nodes_Elem=zeros(nElements,4);
stressXY_nodes_Elem=zeros(nElements,4);

strainXX_nodes_Elem=zeros(nElements,4);
strainYY_nodes_Elem=zeros(nElements,4);
strainXY_nodes_Elem=zeros(nElements,4);

ngaussp=ngp;
Node_Stress=zeros(length(Coord(:,1)),3);
Node_Strain=zeros(length(Coord(:,1)),3);
Node_NoElem=zeros(length(Coord(:,1)),1);
for i=1:nElements
    ue=Eds(i,:);
    eps = [ typeAnalysis E v ];

    [stress, strain, stress_nodes, strain_nodes] = plan4bilinStressStrain...
        (ex(i,:),ey(i,:),eps,ue,ngaussp);

    stressXX_nodes_Elem(i,:)=stress_nodes(1,:);
    stressYY_nodes_Elem(i,:)=stress_nodes(2,:);
    stressXY_nodes_Elem(i,:)=stress_nodes(3,:);
    
    strainXX_nodes_Elem(i,:)=strain_nodes(1,:);
    strainYY_nodes_Elem(i,:)=strain_nodes(2,:);
    strainXY_nodes_Elem(i,:)=strain_nodes(3,:);
    
    % Sum stresses on each node of each element
    for j=1:4
        % Stresses XX
        Node_Stress(Edof(i,1+2*j)/2,1)=Node_Stress(Edof(i,1+2*j)/2,1)+...
                                       stressXX_nodes_Elem(i,j);
        % Stresses YY
        Node_Stress(Edof(i,1+2*j)/2,2)=Node_Stress(Edof(i,1+2*j)/2,2)+...
                                       stressYY_nodes_Elem(i,j);
        % Stresses XY
        Node_Stress(Edof(i,1+2*j)/2,3)=Node_Stress(Edof(i,1+2*j)/2,3)+...
                                        stressXY_nodes_Elem(i,j);
        
        % Strains XX
        Node_Strain(Edof(i,1+2*j)/2,1)=Node_Strain(Edof(i,1+2*j)/2,1)+...
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
    for j=1:4
        node=Edof(i,1+2*j)/2;
        stressXX_nodes_Elem(i,j)=Node_Stress(node,1);
        stressYY_nodes_Elem(i,j)=Node_Stress(node,2);
        stressXY_nodes_Elem(i,j)=Node_Stress(node,3);
        
        strainXX_nodes_Elem(i,j)=Node_Strain(node,1);
        strainYY_nodes_Elem(i,j)=Node_Strain(node,2);
        strainXY_nodes_Elem(i,j)=Node_Strain(node,3);
    end
end
