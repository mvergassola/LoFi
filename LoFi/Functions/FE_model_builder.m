function [Eigenfrequencies_outputs] = FE_model_builder(savefile,gamma)
%   [ Nodes,Elements,Modes,Eigenfrequencies] = FE_model_builder( savefile )
%   FE_model_builder : Executes all the scripts related to the finite element
%   model
%
%   SAVEFILE : matrix of all the inputs
%   Modes : Matrix of the eigenmodes
%   Eigenfrequencies : Vector of the eigenfrequencies -  Hz
%   K_glob : Global Stiffness Matrix
%   M_glob : Global Mass Matrix
%   C_rayleigh : Global Rayleigh damping Matrix
%   C_drag : Global drag damping Matrix
load (savefile)
%--------------------------------------------------------------------------
%Creation of the geometry
disp('Nodes and elements definition : in progress')
[InputsLFM] = geomLFM(savefile);
[Nodes,Elements,~,~] = FE_model_beam(savefile,InputsLFM);
%--------------------------------------------------------------------------
%Creation of additional nodes and elements to increase accuracy.
for i = 1:size(Elements,1)
    Elements(i,1) = i;
end
for i = 1:size(Elements,1)
    [Nodes,Elements] = accuracy_nodes(L_max,Nodes,Elements,i);
end
for i = 1:size(Elements,1)
    Elements(i,1) = i;
end
disp('Nodes and elements definition : completed')
%--------------------------------------------------------------------------
%Creation of the matrices
disp('Creation of the stiffness and mass matrices : in progress')
[K_tot,M_tot,~,~,~] = matrix_assemble(savefile,Nodes,Elements,InputsLFM,gamma);
disp('Creation of the stiffness and mass matrices : completed')
%--------------------------------------------------------------------------
% Soil boundaries conditions definition
[~,K_glob,M_glob] = soil_boundary_condition(K_tot,M_tot);
%--------------------------------------------------------------------------
% Modal analysis
[Modes,Eigenfrequencies] = eigenvalue_calculator(K_glob,M_glob);
disp('Modal mass and stiffness calculation : in progress')
M_mod = diag(Modes'*M_glob*Modes);
K_mod = diag(Modes'*K_glob*Modes);
Eigenfrequencies_outputs(1,:) = { 'Modes' 'Eigenfrequencies - Hz' 'Modal Mass - Kg' 'Modal Stiffness - Nm2'};
for i = 1 : size(K_mod,1)
    Eigenfrequencies_outputs(i+1,:) = {i Eigenfrequencies(i) M_mod(i) K_mod(i)};
end
disp('Modal mass and stiffness calculation : completed')    
disp('Analysis completed')
end