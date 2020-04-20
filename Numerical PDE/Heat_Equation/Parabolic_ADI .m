

clear;
format long;

L = 1; % The length 
W = 1; % The width 
N = 30;% Number of the steps in the x direction 
M = 40;% Number of the steps in the y direction 
alpha = 0.1; % The conductivity 
Tmax = 3.0;
Kmax = 40;

x = linspace(0,L,N);
y = linspace(0,W,M);

dx = L/(N-1); % The step size in the x direction 
dy = W/(M-1); % The step size in the y direction 

dt = Tmax/(Kmax-1); % The step size in the time 

d1 = 0.5*alpha*dt/(dx^2);
d2 = 0.5*alpha*dt/(dy^2);

U = ones(N,M); %initial condition
U(1:N,1) = 0; %bc
U(1:N,M) = 0; %bc
U(1,1:M) = 0; %bc
U(N,1:M) = 0; %bc

%V is for after "stage 1", U->V
V = U; 
%Z is for after "stage 2"  V->Z
Z = U;

%Diagonals to be passed to the FDE matrix
main_diag_x = (1+2*d1)*ones(N-2,1);
off_diag_x = -d1*ones(N-2,1);
main_diag_y = (1+2*d2)*ones(M-2,1);
off_diag_y = -d2*ones(M-2,1);

for k = 1:(Kmax-1)
    for j = 2:(M-1)
        for i=2:(N-1)
            vx(i-1) = d2*U(i,j+1) + (1-2*d2)*U(i,j) + d2*U(i,j-1);
        end
        
        vx(1) = vx(1) + d1*U(1,j);
        vx(N-2) = vx(N-2) + d1*U(N,j);
        V(2:N-1, j) = trid(N-2,main_diag_x,off_diag_x, off_diag_x,vx);
    end
    
    for i=2:(N-1)
        for j=(2:M-1)
           vy(j-1) = d1*V(i+1,j) + (1-2*d1)*V(i,j) +d1*V(i-1,j); 
        end
        vy(1) = vy(1) + d2*V(i,1);
        vy(M-2) = vy(M-2) + d2*V(i,M);
        Z(i,2:M-1) = trid(M-2, main_diag_y, off_diag_y, off_diag_y, vy);
    end
    
    U = Z;
end

surf(Z);