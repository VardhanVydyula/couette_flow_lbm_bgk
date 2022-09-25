%% Parameters
clc 
clear all
close all

Nx = 3; % x- discritisation
Ny = 9; % y - discritisation
Nc = 9; %D2Q9

cx = [0 1 0 -1 0 1 -1 -1 1]; % x-components of unit vectors along 9 directions 
cy = [0 0 1 0 -1 1 1 -1 -1]; % y-components of unit vectors along 9 directions

w = [4/9 1/9 1/9 1/9 1/9 1/36 1/36 1/36 1/36]; % weights of each direction

dt = 1; % non-dimensional time step
dx = 1; % non-dimensional space step
c = 1; % non-dimensional velocity step = dx/dt
tau = 0.9*dt; % relaxation time
u_w = 0.1*c; % Top wall velocity
L2 = 1; % L2 error
count = 1;
f = ones(Nx,Ny,Nc); 
%% Initialisation

for i = 1:1:length(w)
    f(:,:,i) = w(i)*1; % f_eq(rho = 1, u = [0,0] for all (x,y))
end

%% Solution

[rho, u] = macro_prop(f, cx, cy);

while (L2 > 10^(-4) && count < 5000)
    f_new = couette_lbmbgk(f, rho, u, w, dt, tau, u_w, cx, cy, c);
    [rho, u] = macro_prop(f_new, cx, cy);
    u_analy = (u_w/Ny)*(0.5:1:(Ny-0.5));
    %u_mag(:,:) = sqrt(u(:,:,1).^2 + u(:,:,2).^2); 
    L2 = sum((u(2, :, 1) - u_analy).^2);
    f = f_new;
    count =count+ 1;
end

%%
plot(u(2,:,1), 0.5:1:(Ny-0.5), 'ro', 'LineWidth', 2)
hold on ;
%plot(u(2,:,1), 0.5:1:(Ny-0.5), 'ko')
%plot(u(3,:,1), 0.5:1:(Ny-0.5), 'g+')
%close all
%plot (linspace(1,50,50), u_analy)

plot(u_analy, 0.5:1:(Ny-0.5), 'b-', 'LineWidth', 1)
legend('u_{LBM}', 'u_{analytical}', 'FontSize', 18)
xlabel('Veloctiy (u)', 'FontSize', 18);
ylabel('Channel Height (y)', 'FontSize', 18)
title('Comparison of analytical and LBM simulated velocity profile in couette flow ', 'FontSize', 18)
%plot (linspace(1,50,50), u(2,:, 1), 'ro')