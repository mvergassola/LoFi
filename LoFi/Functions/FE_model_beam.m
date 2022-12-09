function [Nodes,Elements,top_node,Nodes_fixed] = FE_model_beam(savefile,InputsLFM)
%   [Nodes,Elements,Nodes_per_legs] = FE_model_monopile(savefile)
%   FE_model_monopile : Creates a finite element model of the substructure 
%   for a monopile based foundation.
%   Nodes : [ Nodes ID , x , y , z] 
%   Elements : [ElID , Node 1 , Node 2 , D1 , D2 , t1 , t2 , E , rho , nu , Ca , Cd , marine growth , Type of element (1 = leg, 2 = brace , 3 = tower, 4 = pile sleeve, 5 = TP, 6 = rigid link, 7 = flexible joint, 8=foundation pile(s), 9=monopile), Material ID] 
%   SAVEFILE : matrix of all the inputs
%   InputsLFM : matrix of the geometric inpunts of the LFM

load (savefile,'E_s','rho_s','nu_s','Nb','Jh')
%==========================================================================
% Definition of the first model node
cpt_nodes = 1;
Nodes_x = 0;
Nodes_y = 0;
Nodes(1,:) = [cpt_nodes, Nodes_x, Nodes_y, 0];
cpt_nodes = cpt_nodes + 1;
% Definition of the geometry of the hollow beam
Monopile_nb_nodes = Nb+1;
L_elem_Monopile = (Jh)/(Monopile_nb_nodes - 1);
origin_diameter = InputsLFM(1)*2;
origin_thickness = InputsLFM(2);
% Definition of the other model nodes
z = L_elem_Monopile;
cpt_elem = 1;
for i=1:Monopile_nb_nodes-1
   Nodes(cpt_nodes,:) = [cpt_nodes , Nodes_x , Nodes_y , z];
   Elements(cpt_elem,1:10) = [cpt_elem, cpt_nodes-1, cpt_nodes, origin_diameter, origin_diameter, origin_thickness, origin_thickness, E_s, rho_s, nu_s];
   Elements(cpt_elem,14) = 9;
   Elements(cpt_elem,15) = 1;   
   cpt_elem = cpt_elem+1;
   z = z + L_elem_Monopile;
   cpt_nodes = cpt_nodes+1;
end
% Definition of the top and fixed (bottom) nodes
top_node = cpt_nodes-1;
Nodes_fixed(1) = 1;
Nodes_fixed(2:4) = 0;
end