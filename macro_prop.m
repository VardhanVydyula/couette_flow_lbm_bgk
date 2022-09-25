function [rho, u] = macro_prop(f, cx, cy)

[Nx,Ny,Nc] = size(f);

rho = zeros(Nx, Ny);
u = zeros(Nx, Ny, 2);

rho(:,:) = sum(f,3);

for i = 1:1:Nx
    for j = 1:1:Ny
        for k = 1:1:Nc
            u(i,j,1) = u(i,j,1) + (f(i,j,k)*(cx(k)))/rho(i,j);
            u(i,j,2) = u(i,j,2) + (f(i,j,k)*(cy(k)))/rho(i,j);
        end
    end
end
