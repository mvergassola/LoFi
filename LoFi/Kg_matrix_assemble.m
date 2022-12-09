function [K_geo_tot,K_elem_loc_export] = Kg_matrix_assemble(Nodes,Elements,disp_elem_loc)
%  [ K_tot , M_tot ] = matrix_assemble(Nodes, Elements )
%   matrix_assemble : This function transforms the elements of the geometric
%   stiffness matrice from the local frame of reference to the global frame
%   of reference. Then the total geometric stiffness matrix is defined.
%
%   K_tot : the total stiffness matrix of the system
%   disp_elem_loc: vector containing the nodal displacements in the local
%   frame of reference.
%   Nodes : [ Nodes ID , x , y , z] 
%   Elements : [ElID , Node 1 , Node 2 , D1 , D2 , t1 , t2 , E , rho , nu , Ca , Cd , marine growth , Type of element (1 = leg, 2 = brace , 3 = tower), Material ID , Flooded (1=Yes 0=No)] 
%   SAVEFILE : matrix of all the inputs

K_geo_tot = zeros(6*size(Nodes,1));
index = 0;
% The local frame of reference is explained in the global frame of reference

Z= [0 0 1];
Y= [0 1 0];
% figure
for k=1:size(Elements,1)
        Node1= Elements(k,2);
        Node2= Elements(k,3);     
    
    if Nodes(Node1,2) == Nodes(Node2,2) && Nodes(Node1,3) == Nodes(Node2,3)
        x1 = (Nodes(Node2,2:4) - Nodes(Node1,2:4)); 
        y1 = cross(x1,Y);
        z1 = cross(x1,y1);
    else   
        x1 = (Nodes(Node2,2:4) - Nodes(Node1,2:4)); 
        y1 = cross(x1,Z);
        z1 = cross(x1,y1);
    end
        x1 = x1/norm(x1);    
        y1 = y1/norm(y1);
        z1 = z1/norm(z1); 
        
%Transition matrix

T = [x1;y1;z1]; % The transfer martix 3x3
T = blkdiag(T,T,T,T); % The 12x12 transfer matrix

% The elements matrices in the global frame of reference
[K_geo_elem_loc,force_multiplier,L] = geometric_stiffness(Nodes,Elements,k);
K_geo_elem_loc = force_multiplier/L*(disp_elem_loc((k-1)*6+1,1,2)-disp_elem_loc((k-1)*6+1,1,1))*K_geo_elem_loc;
K_geo_elem_glob = T'*K_geo_elem_loc*T;

K_geo_elem_glob = (K_geo_elem_glob + K_geo_elem_glob')/2;

K_geo_tot(1+6*(Node1-1):6*Node1,1+6*(Node1-1):6*Node1) = K_geo_elem_glob(1:6,1:6)   + K_geo_tot(1+6*(Node1-1):6*Node1,1+6*(Node1-1):6*Node1);
K_geo_tot(1+6*(Node1-1):6*Node1,1+6*(Node2-1):6*Node2) = K_geo_elem_glob(1:6,7:12)  + K_geo_tot(1+6*(Node1-1):6*Node1,1+6*(Node2-1):6*Node2);
K_geo_tot(1+6*(Node2-1):6*Node2,1+6*(Node1-1):6*Node1) = K_geo_elem_glob(7:12,1:6)  + K_geo_tot(1+6*(Node2-1):6*Node2,1+6*(Node1-1):6*Node1);
K_geo_tot(1+6*(Node2-1):6*Node2,1+6*(Node2-1):6*Node2) = K_geo_elem_glob(7:12,7:12) + K_geo_tot(1+6*(Node2-1):6*Node2,1+6*(Node2-1):6*Node2);

K_elem_loc_export(index*6+1:(index+2)*6,index*6+1:(index+2)*6) = K_geo_elem_loc;
index = index + 2;
end