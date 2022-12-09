              %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
              %                 LoFi                   %
              %  Low-fidelity model of a full lattice  %
              %     wind turbine support structure     %
              %            Marco Vergassola            %
              %              December 2022             %
              %                V 1.0.0                 %
              %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars; close all;

%%%%%%%%%%%%%%%%%%%%%%%%% INPUTS DEFINITION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Enter the configuration name :
Inputs = 'Inputs.mat' ;

%==========================================================================
%==========================================================================
% Section 1
% Material Properties Inputs 

rho_s= 7850;       % Steel density - kg/m3
E_s = 2.1e11;    % Steel Young's modulus - Pa
nu_s = 0.3;        % Steel Poisson's ratio

% Provide the material properties for the design of the element connecting
% the structure top and the RNA.

rho_RNA = 0.785;   % Connecting element density - kg/m3
E_RNA = 2.1e14;    % Connecting element Young's modulus - Pa
nu_RNA = 0.3;      % Connecting element Poisson's ratio

%==========================================================================
%==========================================================================
% Section 2
% Structure Geometry Inputs

L_max = 100.0;      % Maximum length of an element in the FE model - meter

Jh = 180;      % Lattice structure Height - meter
Nb = 10;        % Number of bay
L_bottom = 18;      % Bottom width - meter
L_top = 18;        % Top width - meter

D_leg = 2;  % Legs diameter - meter (outer diameter)
t_leg = 2/50; % Legs wall thickness - meter

%--------------------------------------------------------------------------
% Section 2.1
% Bracing Inputs
Brace_pattern = 'X';     % The bracing type :
                                        % 'X' : X bracing
                                        % '?' : other bracings to be added

D_brace = 1;         % Brace diameter - meter
t_brace = 1/40;         % Brace thickness - meter

%==========================================================================
%==========================================================================
% Section 3
% Soil stiffness definition

% Stiffness_type : 
             % 0 : Rigid foundation
             % 1 : Other options coming soon

stiffness_type = 0;

%%%%%%%%%%%%%%%%%%%%%% END OF THE INPUTS DEFINITION %%%%%%%%%%%%%%%%%%%%%%%
%% FE model
% !!! DO NOT MODIFY THIS SECTION !!!
save (Inputs,'L_max','Nb','Jh','L_bottom','L_top','D_leg','t_leg',...
    'D_brace','t_brace','rho_s','E_s','nu_s','rho_RNA','E_RNA','nu_RNA',...
    'stiffness_type');

[gamma] = gammaPredictor(L_bottom,Jh/Nb,Nb,D_leg,D_brace);
[Eigenfrequencies_outputs] = FE_model_builder(Inputs,gamma);