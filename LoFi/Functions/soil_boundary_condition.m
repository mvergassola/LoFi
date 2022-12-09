function [K_foundation,K_glob,M_glob] = soil_boundary_condition(K_tot,M_tot)
%   [K_foundation,K_glob,M_glob,C_drag] = soil_boundary_condition(savefile,Nodes,Elements,K_tot,M_tot,C_drag,Weight_marine_growth_tot,Boyancy_force,Nodes_fixed,LJF)
%   soil_boundary_condition : Define soil boundary condition for the model 
%
%   K_tot : The total stiffness matrix
%   M_tot : The total mass matrix
%   K_glob : the global stiffness matrix
%   M_glob : the global mass matrix

disp('Soil boundaries conditions definition : in progress')

% Definition of the fixed nodes
%--------------------------------------------------------------------------
% Rigid foundation
K_glob = K_tot(7:end, 7:end); 
M_glob = M_tot(7:end, 7:end);
K_foundation = 0;
disp('Soil boundaries conditions definition : completed')
end