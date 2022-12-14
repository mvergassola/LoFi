function [K_tot,M_tot,K_elem_glob_export,M_elem_glob_export,T_export] = matrix_assemble(savefile,Nodes,Elements,InputsLFM,gamma)
%  [ K_tot , M_tot ] = matrix_assemble(Nodes, Elements )
%   matrix_assemble : This function transforms the elements matrices from 
%   the local frame of reference to the global frame of reference.  
%   Then the total mass matrix and the total stiffness matrix are defined.
%
%   K_tot : the total stiffness matrix of the system
%   M_tot : the total mass matrix of the system
%   Nodes : [ Nodes ID , x , y , z] 
%   Elements : [ElID , Node 1 , Node 2 , D1 , D2 , t1 , t2 , E , rho , nu , Ca , Cd , marine growth , Type of element (1 = leg, 2 = brace , 3 = tower), Material ID , Flooded (1=Yes 0=No)] 
%   SAVEFILE : matrix of all the inputs

%Initialize stiffness and mass matrices
K_tot = zeros(6*size(Nodes,1));
M_tot = zeros(6*size(Nodes,1));
index = 0;
% The local frame of reference is explained in the global frame of reference
Z= [0 0 -1];
Y= [1 0 0];

for k=1:size(Elements,1 )
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
    [K_elem_loc, M_elem_loc] = element_matrices(savefile,Nodes,Elements,k,gamma,InputsLFM);        

    M_elem_glob = T'*M_elem_loc*T;
    K_elem_glob = T'*K_elem_loc*T;

    K_elem_glob = (K_elem_glob + K_elem_glob')/2;
    M_elem_glob = (M_elem_glob + M_elem_glob')/2;

    K_tot(1+6*(Node1-1):6*Node1,1+6*(Node1-1):6*Node1) = K_elem_glob(1:6,1:6)   + K_tot(1+6*(Node1-1):6*Node1,1+6*(Node1-1):6*Node1);
    K_tot(1+6*(Node1-1):6*Node1,1+6*(Node2-1):6*Node2) = K_elem_glob(1:6,7:12)  + K_tot(1+6*(Node1-1):6*Node1,1+6*(Node2-1):6*Node2);
    K_tot(1+6*(Node2-1):6*Node2,1+6*(Node1-1):6*Node1) = K_elem_glob(7:12,1:6)  + K_tot(1+6*(Node2-1):6*Node2,1+6*(Node1-1):6*Node1);
    K_tot(1+6*(Node2-1):6*Node2,1+6*(Node2-1):6*Node2) = K_elem_glob(7:12,7:12) + K_tot(1+6*(Node2-1):6*Node2,1+6*(Node2-1):6*Node2);

    M_tot(1+6*(Node1-1):6*Node1,1+6*(Node1-1):6*Node1) = M_elem_glob(1:6,1:6)   + M_tot(1+6*(Node1-1):6*Node1,1+6*(Node1-1):6*Node1);
    M_tot(1+6*(Node1-1):6*Node1,1+6*(Node2-1):6*Node2) = M_elem_glob(1:6,7:12)  + M_tot(1+6*(Node1-1):6*Node1,1+6*(Node2-1):6*Node2);
    M_tot(1+6*(Node2-1):6*Node2,1+6*(Node1-1):6*Node1) = M_elem_glob(7:12,1:6)  + M_tot(1+6*(Node2-1):6*Node2,1+6*(Node1-1):6*Node1);
    M_tot(1+6*(Node2-1):6*Node2,1+6*(Node2-1):6*Node2) = M_elem_glob(7:12,7:12) + M_tot(1+6*(Node2-1):6*Node2,1+6*(Node2-1):6*Node2);

    K_elem_glob_export(index*6+1:(index+2)*6,index*6+1:(index+2)*6) = K_elem_loc;
    M_elem_glob_export(index*6+1:(index+2)*6,index*6+1:(index+2)*6) = M_elem_loc;
    T_export(index*6+1:(index+2)*6,index*6+1:(index+2)*6) = T;
    index = index + 2;
end

end
