% Heat_Equation_Explicit.m

% Preparation of the Working Station 
clc;
close all,
clear 

% The 1D Problem Specification 
L=1.0; % The length of the wire 
T=1.0; % The Time durition when the heat is conducted through the length 

% Creating the 1D mesh that each node represents a solution 
nX= 100; % The number of the nodes 
dx=L/nX; % The length of each spacial step 

nT=2000; % The number of the time steps 
dt=T/nT; % The length of each time step 

con=1.5; % The wire conductivity coefficient 
alpha=con*(dt/dx*dx);

% Initial Value Parameters 
for i=1:nX+1
    x(i)=(i-1)*dx;
    u(i,1)=sin(pi*x(i));
end

% The Boundary Condition where T=0
for n=1:nT+1
   u(1,n)=0;
   u(nX+1,n)=0;
   time(n)=(nT-1)*dt;
end

% Implementation of the Explicit Method 
 for n=1:nT
     for i=2:nX
         u(i,n+1)=u(i,n)+ alpha*(u(i+1,n)-2*u(i,n)+u(i-1,n));
     end
 end
 
 % Graphical representation of the temperature at different selected times
figure(1)
plot(x,u(:,1),'-',x,u(:,100),'-',x,u(:,300),'-',x,u(:,600),'-')
title('Temperature within the explicit method')
xlabel('X')
ylabel('T')
figure(2)
mesh(x,time,u')
title('Temperature within the explicit method')
xlabel('X')
ylabel('Temperature')
      


