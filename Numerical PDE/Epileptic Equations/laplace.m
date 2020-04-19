% laplace.m

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
% Control by plotting the mesh 
plot(x,y);
beta=dx/dy;

% Specification of the method applied 
% The Iterative Methods are symbolized by the specification of " Tolerance"
Tolerance=1e-5;
converged=0; % The convergence of the iteration is set to be false
n_Iter=0; % The number of the iterations are stored 
U=ones(N,N); % The solution matrix 
V=U; % The updated matrix after the iteration is represented by V , storing the most updated ones

% The computational Iteration on the mesh 
while (~converged)
    n_Iter = n_Iter + 1;
    for j = 2:N-1
        for i=2:N-1
           term1 = U(i+1,j)+U(i-1,j);
           term2 = U(i,j+1)+U(i,j-1);
           V(i,j) = 0.5*(term1+beta^2*term2)/(1+beta^2);
        end
    end
     
    error = max(max(abs(V-U)));
    if (error < Tolerance)
        converged = 1;
    end
    
    % increment solution
    U=V;
end

surf(x,y,V);