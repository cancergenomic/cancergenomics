%Laplace_SOR.m

% Preparing the workspace 
clear;
close all;
format long;

% Now , construct the 2D mesh for the "stencil" to move around 
% The " Unit Square" is defined as:

L = 1.0; % The length of the mesh which represents the x-axis 
W = 1.0; % The width of the mesh which represents the y-axis 

N=100; % The number of the nodes in each direction x , y 
% The x-axis 
x=linspace(1,L,N);
% The y-axis 
y=linspace(1,W,N);
% The diffrential of the space 
dx=L/(N-1);
dy=L/(N-1);
beta=dx/dy;
% Control by plotting the mesh 
%plot(x,y);

% Make some specification for the iterative method 
omega=1.9;
tolerance=1e-4;
U=ones(N,N);
V=U;
W=U;

converged=0;
n_Iter=0;

while(~converged)
    n_Iter=n_Iter+1;
    for j=2:N-1
        for i=2:N-1
            term1=U(i+1,j)+V(i-1,j);
            term2=U(i,j+1)+V(i,j-1);
            V(i,j) = 0.5*(term1+beta^2*term2)/(1+beta^2);
            W(i,j) = (1-omega)*U(i,j)+omega*V(i,j);
        end
    end
   error= max(max(abs(V-U)));
   if error < tolerance
       converged=1;
   end
   U=W;
end

surf(x,y,V);