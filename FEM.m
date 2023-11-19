function [U,U_exact] = FEM(N)

x1 = linspace(0, 1, N+1);
y1 = linspace(0, 1, N+1);
x = linspace(0, 1, N);
y = linspace(0, 1, N);
[X, Y] = meshgrid(x1, y1);
h  = 1/ N;
h1 = 1/(N+1);
A = sparse(zeros((N+1)^2,(N+1)^2));
U = sparse(zeros((N+1)^2,1));
K = sparse(zeros((N+1)^2,1));
local_k = zeros(4,1);
local_a = zeros(4,4);

% Define the source term and the exact solution
f = @(x, y) 400.*(x.^4-y.^4).*sin(20.*x.*y);
U_exact = (X.^2 - Y.^2) .* sin(20.*X.*Y);

% Define u_x and u_y of the boundary condition
g1 = @(x) (20.*x.*(x.^2-1).*cos(20.*x)-2.*sin(20.*x));
g2 = @(y) (20.*y.*(1-y.^2).*cos(20.*y)+2.*sin(20.*y));

% Define the basis function and their gradiant

phi{1} = @(x,y) (1 - x) .* (1 - y);
phi{2} = @(x,y) x .* (1 - y);
phi{3} = @(x,y) x .* y;
phi{4} = @(x,y) (1 - x) .* y;

% phi{1} = @(x, xi, y, yj) (1 - (x-xi)./h) .* (1 - (yj-y)./h);
% phi{2} = @(x, xi, y, yj) (x-xi)./h .* (1 - (y-yj)./h);
% phi{3} = @(x, xi, y, yj) (xi-x)./h .* (yj-y)./h;
% phi{4} = @(x, xi, y, yj) (1 - (xi-x)./h) .* (yj-y)./h;

% Define the shape function derivatives
dN1_dx = @(x, y) -1 + y;
dN1_dy = @(x, y) -1 + x;
dN2_dx = @(x, y) 1 - y;
dN2_dy = @(x, y) -x;
dN3_dx = @(x, y) y;
dN3_dy = @(x, y) x;
dN4_dx = @(x, y) -y;
dN4_dy = @(x, y) 1-x;

% dN1_dx = @(x, xi, y, yj) (-1 + (yj-y)./h)./h;
% dN1_dy = @(x, xi, y, yj) (1 - (x-xi)./h)./h;
% 
% dN2_dx = @(x, xi, y, yj) (1 - (y-yj)./h)./h;
% dN2_dy = @(x, xi, y, yj) -(x-xi)./h./h;
% 
% dN3_dx = @(x, xi, y, yj) (y-yj)./h./h;
% dN3_dy = @(x, xi, y, yj) (x-xi)./h./h;
% 
% dN4_dx = @(x, xi, y, yj) (yj-y)./h./h;
% dN4_dy = @(x, xi, y, yj) -(1 - (xi-x)./h)./h;

% Define function handles for derivatives
dN_dx = {dN1_dx, dN2_dx, dN3_dx, dN4_dx};
dN_dy = {dN1_dy, dN2_dy, dN3_dy, dN4_dy};

for i = 1:N
    for j = 1:N
        global_indices  = [i*(N+1) + j, i*(N+1) + j+1, (i-1)*(N+1) + j+1,(i-1)*(N+1) + j];
        for p = 1:4
            integrand2 = @(x,y) phi{p}(x, y) .* f(x,y);
%             integrand2 = @(x,y) phi{p}(x, x(i), y, y(j)) .* f(x,y);
%             local_k(p) = integral2(integrand2,(i-1)*h,i*h,(j-1)*h,(j)*h);
            %local_k(p) = integral2(integrand2,0,1,0,1);
            local_k(p) = integrand2((i-0.5)*h,(j-0.5)*h)*h*h;
            k = global_indices(p);

            if i > 1 && j>1 && i < N && j<N
               K(k) = K(k) + local_k(p);
            else
                if i == 1 || j == 1
                    K(k) = 0;
                    A(k,:) = 0;
                    A(:,k) = 0;
                    A(k,k) = 1;
                elseif j == N && i < N
                    integrand3 = @(x) phi{p}(x,1) .* g1(x);
%                   integrand3 = @(x) phi{p}(x, x(i), 1, y(j)) .* g1(x);
%                   K(k) = K(k) +  local_k(p) + integral(integrand3, (i-1)*h,(i)*h);
                    %K(k) = K(k) +  local_k(p) + integral(integrand3,0,1);
                    K(k) = K(k) + local_k(p) + integrand3((i-0.5)*h)*h;
                elseif i == N 
                    integrand4 = @(y) phi{p}(1,y) .* g2(y);
%                   integrand4 = @(y) phi{p}(1, x(i), y, y(j)) .* g2(y);
%                   K(k) = K(k) + local_k(p) + integral(integrand4,(j-1)*h,(j)*h);
                    %K(k) = K(k) + local_k(p) + integral(integrand4,0,1);
                    K(k) = K(k) + local_k(p) + integrand4((j-0.5)*h)*h;
                end
            end

            for q = 1:4
                integrand1 = @(x, y) dN_dx{p}(x, y) .* dN_dx{q}(x, y) + ...
                    dN_dy{p}(x, y) .* dN_dy{q}(x, y);
%                 integrand1 = @(x, y) dN_dx{p}(x, xi(p), y, yj(p)) .* dN_dx{q}(x, xi(q), y, yj(q)) + ...
%                     dN_dy{p}(x, xi(p), y, yj(p)) .* dN_dy{q}(x, xi(q), y, yj(q));
                local_a(p,q) = integrand1((i-0.5)*h,(j-0.5)*h);
                %local_a(p,q) = integral2(integrand1, 0, 1,0, 1);
%                 local_a(p,q) =  integral2(integrand1,(i-1)*h,i*h,(j-1)*h,(j)*h)/(h^2);
                n = global_indices(q);
                if i > 1 && j > 1
%                     A(m,n) = A(m,n) + local_a(p,q);
                   A(k,n) = A(k,n) + local_a(p,q);
                end
            end
        end
    end
end

% for j = 1:N
%     for i = 1:N
%         global_indices = [(i-1)*(N+1) + j, (i-1)*(N+1) + j+1, i*(N+1) + j+1, i*(N+1) + j];
%         for p = 1:4
%             integrand2 = @(x, y) phi{p}(x,y) .* f(x,y);
%             %local_k(p) = integral2(integrand2, (i-1)*h, (i)*h,(j-1)*h, (j)*h);
%             local_k(p) = integrand2((i-0.5)*h,(j-0.5)*h)*h*h;
%             k = global_indices(p);
%             if i > 1 && j>1 && i < N && j<N
%                K(k) = K(k) + local_k(p);
%             elseif i == 1 || j == 1
%                 K(k) = 0;
%                 A(:,k) = 0;
%                 A(k,:) = 0;
%                 A(k,k) = 1;
%             elseif j == N
%                 integrand3 = @(x) phi{p}(x,1) .* g1(x);
%                 %K(k) = K(k) +  local_k(p) + integral(integrand3, 0, 1);
%                 K(k) = K(k) +  local_k(p) + integrand3((i-0.5)*h)*h;
%             elseif i == N
%                 integrand4 = @(y) phi{p}(1,y) .* g2(y);
%                 %K(k) = K(k) + local_k(p) + integral(integrand4, 0, 1);
%                 K(k) = K(k) + local_k(p) + integrand4((j-0.5)*h)*h;
%             end
%         end
%     end
% end

U = full(A \ K);
U = reshape(U,N+1,N+1);
end




