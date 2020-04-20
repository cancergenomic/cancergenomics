% Heat_Equation_Implicit.m

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

% Implementatin of the Implicit Methods 

aa(1:nX-2)=-alpha; % Off diagonal Elements 
bb(1:nX-1)=2*alpha+1; % The main diagonal elements 
cc(1:nX-2)=-alpha; % The Off Diagonal elements 

MM=inv(diag(bb,0)+diag(aa,-1)+diag(cc,+1));

for n=2:nT
    uu=u(2:nX,nT-1);
    u(2:nX,nT)=MM*uu;
end

% Graphical representation of the temperature at different selected times
figure(1)
plot(x,u(:,1),'-',x,u(:,100),'-',x,u(:,300),'-',x,u(:,600),'-')
title('Temperature within the fully implicit method')
xlabel('X')
ylabel('T')
figure(2)
mesh(x,time,u')
title('Temperature within the fully implicit method')
xlabel('X')
ylabel('Temperature')
