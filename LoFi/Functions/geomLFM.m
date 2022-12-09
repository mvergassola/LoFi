function [InputsLFM] = geomLFM(savefile)
%   geomLFM : compute the geometric inpunts for the low-fidelity model via
%   mechanical and mass homogenization of the space-frame bay
%
%   InputsLFM : matrix of the geometric inpunts of the LFM
%   SAVEFILE : matrix of all the inputs

disp('Low-fidelity model elements definition : in progress')
load (savefile)
%Definition of the parameters used to compute the homogenized geometry
rho = rho_s; %Material density - [kg/m3]
Hb = Jh/Nb; %Height of the bay - [m]
R_o_l = D_leg/2;%External radius hollow section legs - [m] 
R_i_l = R_o_l - t_leg; %Interal radius hollow section legs - [m]
R_o_b = D_brace/2; %External radius hollow section braces - [m] 
R_i_b = R_o_b - t_brace; %Interal radius hollow section braces - [m]
Width = L_bottom; %Width of the bay - [m]
L_l = Hb; %Length of the legs - [m]
L_b = sqrt(L_l^2 + Width^2); %Length of the braces - [m]
dist_l = (Width*sqrt(2))/2;
dist_b = Width/2;

%% Solve homogenization problem 
syms r_o r_i
%Compute cross-section area of legs and braces
A_l = pi*(R_o_l^2 - R_i_l^2);
A_b = pi*(R_o_b^2 - R_i_b^2);
%Compute mass of bay
m_l = rho*A_l*L_l;
m_b = rho*A_b*L_b;
m = 4*m_l + 8*m_b;
%Compute second area moment of inertia of bay
J_l = pi/2*(R_o_l^4 - R_i_l^4) + dist_l^2*A_l;
J = 4*J_l;
%Compute mass and second area moment of inertia of equivalent hollow beam
m_eq = rho*Hb*pi()*(r_o^2 - r_i^2);
J_eq = pi/2*(r_o^4 - r_i^4);
%Construct and solve system of equations
eqn1 = J - J_eq == 0;
eqn2 = m - m_eq == 0;
S = solve([eqn1 eqn2],[r_o r_i], 'ReturnConditions', true);
sol1 = vpa(S.r_o);
sol2 = vpa(S.r_i);
%Compute characteristic parameters of the homogenized section
Ro = double(abs(sol1(1)));
Ri = double(abs(sol2(1)));
t = Ro - Ri;
A_new = pi()*(Ro^2 - Ri^2);
Ix_new = pi()/4*(Ro^4 - Ri^4);
Iy_new = Ix_new;
Iz_new = Ix_new + Iy_new;
%Compute G values in all the directions based on Noor et al (1979)
E = 2.1E+11;
AGx = 4*E*A_b*Hb*(dist_b*2)^2/(L_b^3);
JGz = 2*E*A_b*Hb*(dist_b*2)^4/(L_b^3);
Gx = AGx/A_new;
Gz = JGz/Iz_new;
%Create output matrix
InputsLFM = [Ro, t, Gx, Gz];

end