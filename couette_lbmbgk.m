function f_new = couette_lbmbgk(f, rho, u, w, dt, tau, u_w, cx, cy, c)

[Nx, Ny, Nc] = size(f);

f_temp = zeros(Nx, Ny, Nc);

f_new = zeros(Nx, Ny, Nc);
%% Collision
for i = 1:1:Nx
    for j = 1:1:Ny
        for k = 1:1:Nc
            f_eq(i,j,k) = w(k)*rho(i,j)*(1 + 3*(cx(k)*u(i,j,1)+cy(k)*u(i,j,2))...
                + 4.5*(cx(k)*u(i,j,1)+cy(k)*u(i,j,2))^2-1.5*(u(i,j,1)^2 + u(i,j,2)^2));
            f_temp(i,j,k) = f(i,j,k) - (dt/tau)*(f(i,j,k) - f_eq(i,j,k));
        end
    end
end

%% Streaming

for i = 1:1:Nx
    for j = 1:1:Ny
        for k = 1:1:Nc
            
            try %Try-Catch algorithm because some cases may end up with 0 or negative indices of f, which are later taken care by the BC
                f_new(i+cx(k), j + cy(k), k) = f_temp(i,j,k); 
            catch
                continue;
            end
            
            
            %{
            if j == 1 % bottom wall
                if k ==2
                    f_new(i, j, k) = f_temp(i,j,4);
                elseif k==5
                    f_new(i, j, k) = f_temp(i,j,7);
                elseif k == 6
                    f_new(i, j, k) = f_temp(i,j,8);
                else
                    f_new(i,j,k)= 11111111111111111111;
                end
                
            elseif j == Ny % top wall
                if k==4
                    f_new(i, j, k) = f_temp(i,j,2);
                elseif k ==7
                    f_new(i, j, k) = f_temp(i,j,5) - (1/(6*c))*u_w;
                elseif k ==8
                    f_new(i, j, 8) = f_temp(i,j,6) + (1/(6*c))*u_w;
                else
                  f_new(i,j,k)= 11111111111111111111;
                end
                  
            elseif i == 1 % left end
                 f_new(i, j, 1) = f_temp(Nx,j,1);
                 f_new(i, j, 5) = f_temp(Nx,j,5);
                 f_new(i, j, 8) = f_temp(Nx,j,8);
                 
            elseif i == Nx % right end
                f_new(i, j, 3) = f_temp(1,j,3);
                f_new(i, j, 6) = f_temp(1,j,6);
                f_new(i, j, 7) = f_temp(1,j,7);
                
            end
            %}
        end
    end
end
f_new = f_new(1:end-1, 1:end-1, :); % Discarding all the cases with Nx+1, Ny+1 indices which are taken care by the periodic BC

            % "index + 1" in the code because book has c0 through c8, but
            % matlab has c1 through c9 as the velocity vectors
            
            %Bottom Wall
            f_new(:, 1, 3) = f_temp(:,1,5); %f(x, y1, 2) = f*(x,y1,4)
            f_new(:, 1, 6) = f_temp(:,1,8); %f(x, y1, 5) = f*(x,y1,7)
            f_new(:, 1, 7) = f_temp(:,1,9); %f(x, y1, 6) = f*(x,y1,8)
            %Top Wall (Moving Wall)
            f_new(:, Ny, 5) = f_temp(:,Ny,3);%f(x, yN, 4) = f*(x, yN, 2)
            f_new(:, Ny, 8) = f_temp(:,Ny,6) - (1/(6*c))*u_w;%f(x, yN, 7) = f*(x, yN, 5) - (1/6c)u_w
            f_new(:, Ny, 9) = f_temp(:,Ny,7) + (1/(6*c))*u_w;%f(x, yN, 8) = f*(x, yN, 6) + (1/6c)u_w
            % Right to left periodic
            f_new(1, :, 2) = f_temp(Nx,:,2);%f(x1, y, 1) = f*(xN, y, 1)
            f_new(1, :, 6) = f_temp(Nx,:,6);%f(x1, y, 5) = f*(xN, y, 5)
            f_new(1, :, 9) = f_temp(Nx,:,9);%f(x1, y, 8) = f*(xN, y, 8)
            %Left to Right periodic
            f_new(Nx, :, 4) = f_temp(1,:,4);%f(xN, y, 3) = f*(x1, y, 3)
            f_new(Nx, :, 7) = f_temp(1,:,7);%f(xN, y, 6) = f*(x1, y, 6)
            f_new(Nx, :, 8) = f_temp(1,:,8);%f(xN, y, 7) = f*(x1, y, 7)
