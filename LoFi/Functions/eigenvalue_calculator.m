 function [modes,eigenfrequencies] = eigenvalue_calculator(K_glob,M_glob)
%   [modes,eigenfrequencies] = eigenvalue_calculator(K_glob,M_glob,Nodes,Elements,Nodes_fixed,savefile)
%   eigenvalue_calculator :  Calculate the modes and the eigenvalues of the
%   model.
%
%   K_glob : Global Stiffness Matrix
%   M_glob : Global Mass Matrix
%   modes : Matrix of the eigenmodes
%   eigenfrequencies : Vector of the eigenfrequencies -  Hz
%   Nodes : [ Nodes ID , x , y , z] 
%   Elements : [ElID , Node 1 , Node 2 , D1 , D2 , t1 , t2 , E , rho , nu , Ca , Cd , marine growth , Type of element (1 = leg, 2 = brace , 3 = tower, 4 = pile sleeve), Material ID] 
%   SAVEFILE : matrix of all the inputs
%   Nodes_fixed :  The nodes fixed to the seabed
%   TF_Kp : Array representing the the soil stiffness

%Calculation of the eigenvalues and eigenmodes

disp('Modal analysis : in progress')
[modes,eigenfrequencies] = eig(K_glob,M_glob) ;
eigenfrequencies = sqrt((diag(eigenfrequencies)))/2/pi;
[eigenfrequencies,I]=sort(real(eigenfrequencies));
eigenfrequencies = [ I eigenfrequencies ];

disp('Modal analysis : completed')    
eigenfrequencies = eigenfrequencies(:,2);
end

