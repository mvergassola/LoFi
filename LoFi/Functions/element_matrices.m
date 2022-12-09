function [ K_elem, M_elem] = element_matrices( savefile,Nodes,Elements,Element_ID,gamma,InputsLFM)
%   [ K_elem, M_elem ] = Element_Matrices( Nodes,Element)
%   Element_Matrices return the mass matrix and the stiffness matrix of one
%   tubular element
%
%   Nodes : [ Nodes ID , x , y , z] 
%   Elements : [ElID , Node 1 , Node 2 , D1 , D2 , t1 , t2 , E , rho , nu , Ca , Cd , marine growth , Type of element (1 = leg, 2 = brace , 3 = tower, 4 = pile sleeve), Material ID] 
%   Element_ID : The element that is going to be treated
%   gamma : vector containing the 3 gamma coefficientes
%   SAVEFILE : matrix of all the inputs
load (savefile)

%Extract information from element k
Node1  = Elements(Element_ID,2);
Node2  = Elements(Element_ID,3);

R1 = Elements(Element_ID,4)/2;
R2 = Elements(Element_ID,5)/2;
t1 = Elements(Element_ID,6);
t2 = Elements(Element_ID,7);

Diameter  = (Elements(Element_ID,4) +  Elements(Element_ID,5))/2;
Thickness = (Elements(Element_ID,6) +  Elements(Element_ID,7))/2;

Iy = (pi/64)*((Diameter)^4-(Diameter-2*Thickness)^4);
Iz = (pi/64)*((Diameter)^4-(Diameter-2*Thickness)^4);
Iz = Iz - 1e-5;
Ix = (Iy+Iz);

if Element_ID < Nb+1
    E   = Elements(Element_ID,8)*gamma(1);
else
    E   = Elements(Element_ID,8);
end

rho = Elements(Element_ID,9);
nu  = Elements(Element_ID,10);

A = pi*(Diameter/2)^2 - pi*((Diameter-2*Thickness)/2)^2;
L = norm(Nodes(Node2,2:4)-Nodes(Node1,2:4)); 
ks = 3/4;
As = ks*A; % cross sectional area effective in shear


if Element_ID < Nb+1
    Gx = InputsLFM(3)*gamma(2);
    Gz = InputsLFM(4)*gamma(3);
    Gy = Gz;
else
    G_isotropic = E/(2*(nu+1)); %Shear coefficient
    Gx = G_isotropic;
    Gz = G_isotropic;
    Gy = Gz;
end
G = [Gx; Gy; Gz]; %Shear coefficients

Phi_y = 12*E*Iz/(G(2)*As*L^2); %Shear deformation parameter
Phi_z = 12*E*Iy/(G(3)*As*L^2);

outer_volume = (L*pi/3)*(R1^2+R2^2+R1*R2);
inner_volume = (L*pi/3)*((R1-t1)^2+(R2-t2)^2+(R1-t1)*(R2-t2));
Volume_elem = outer_volume-inner_volume;
m  = Volume_elem*rho;

%% Define entries of the element mass and stiffness matrices
% Element Mass Matrix
M_elem = zeros(12,12);

M_elem(1,1) = 1/3;
M_elem(2,2) = 13/35+6*Iz/(5*A*(L^2));
M_elem(3,3) = 13/35+6*Iy/(5*A*(L^2));
M_elem(4,4) = Ix/(3*A);
M_elem(5,5) = (L^2)/105+2*Iy/(15*A);
M_elem(6,6) = (L^2)/105+2*Iz/(15*A);

M_elem(7,7) = M_elem(1,1);
M_elem(8,8) = M_elem(2,2);
M_elem(9,9) = M_elem(3,3);
M_elem(10,10) = M_elem(4,4);
M_elem(11,11) = M_elem(5,5);
M_elem(12,12) = M_elem(6,6);

M_elem(2,6) = 11*L/210+Iz/(10*A*L);
M_elem(6,2) = M_elem(2,6);

M_elem(3,5) = -11*L/210-Iy/(10*A*L);
M_elem(5,3) = M_elem(3,5);

M_elem(9,11) = 11*L/210+Iy/(10*A*L);
M_elem(11,9) = M_elem(9,11);

M_elem(7,1) = 1/6 ;
M_elem(1,7) = M_elem(7,1) ;

M_elem(8,2) = 9/70-6*Iz/(5*A*(L^2));
M_elem(2,8) = M_elem(8,2);

M_elem(8,6) = 13*L/420-Iz/(10*A*L);
M_elem(6,8) = M_elem(8,6);

M_elem(9,3) = 9/70-6*Iy/(5*A*(L^2));
M_elem(3,9) = M_elem(9,3);

M_elem(9,5) = -13*L/420+Iy/(10*A*L);
M_elem(5,9) = M_elem(9,5);

M_elem(10,4) = Ix/(6*A);
M_elem(4,10) = M_elem(10,4);

M_elem(11,3) = 13*L/420-Iy/(10*A*L);
M_elem(3,11) = M_elem(11,3);

M_elem(11,5) = -(L^2)/140-Iy/(30*A);
M_elem(5,11) = M_elem(11,5);

M_elem(12,2) = -13*L/420+Iz/(10*A*L);
M_elem(2,12) = M_elem(12,2);

M_elem(12,6) = -(L^2)/140-Iz/(30*A);
M_elem(6,12) = M_elem(12,6);

M_elem(12,8) = -11*L/210-Iz/(10*A*L);
M_elem(8,12) = M_elem(12,8);

M_elem = m*M_elem;

%%
%Element Stiffness matrix
K_elem = zeros(12,12);

K_elem(1,1) = E*A/L;
K_elem(2,2) = 12*E*Iz/((L^3)*(1+Phi_y));
K_elem(3,3) = 12*E*Iy/((L^3)*(1+Phi_z));
K_elem(4,4) = G(1)*Ix/L;
K_elem(5,5) = (4+Phi_z)*E*Iy/(L*(1+Phi_z));
K_elem(6,6) = (4+Phi_y)*E*Iz/(L*(1+Phi_y));
K_elem(7,7) = K_elem(1,1);
K_elem(8,8) = K_elem(2,2);
K_elem(9,9) = K_elem(3,3);
K_elem(10,10) = K_elem(4,4);
K_elem(11,11) = K_elem(5,5);
K_elem(12,12) = K_elem(6,6);

K_elem(7,1) = -E*A/L;
K_elem(1,7) = K_elem(7,1);

K_elem(5,3) = -6*E*Iy/((L^2)*(1+Phi_z));
K_elem(3,5) = K_elem(5,3);

K_elem(6,2) = 6*E*Iz/((L^2)*(1+Phi_y));
K_elem(2,6) = K_elem(6,2);

K_elem(8,2) = -12*E*Iz/((L^3)*(1+Phi_y));
K_elem(2,8) = K_elem(8,2);

K_elem(8,6) = -6*E*Iz/((L^2)*(1+Phi_y));
K_elem(6,8) = K_elem(8,6);

K_elem(9,3) = -12*E*Iy/((L^3)*(1+Phi_z));
K_elem(3,9) = K_elem(9,3);

K_elem(9,5) = 6*E*Iy/((L^2)*(1+Phi_z));
K_elem(5,9) = K_elem(9,5);

K_elem(10,4) = -G(1)*Ix/L;
K_elem(4,10) = K_elem(10,4) ;

K_elem(11,3) = -6*E*Iy/((L^2)*(1+Phi_z));
K_elem(3,11) = K_elem(11,3);

K_elem(11,5) = (2-Phi_z)*E*Iy/(L*(1+Phi_z));
K_elem(5,11) = K_elem(11,5);

K_elem(11,9) = 6*E*Iy/((L^2)*(1+Phi_z));
K_elem(9,11) = K_elem(11,9);

K_elem(12,2) = 6*E*Iz/((L^2)*(1+Phi_y));
K_elem(2,12) = K_elem(12,2);

K_elem(12,6) = (2-Phi_y)*E*Iz/(L*(1+Phi_y));
K_elem(6,12) = K_elem(12,6);

K_elem(12,8) = -6*E*Iz/((L^2)*(1+Phi_y));
K_elem(8,12) = K_elem(12,8);

end

