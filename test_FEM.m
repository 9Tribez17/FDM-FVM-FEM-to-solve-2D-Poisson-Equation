clc;
clear;
% Define grid sizes for convergence study
grid_sizes = (10:10:100);
errors = [];
count = 0;
%Loop over the grid sizes to compute errors
for n = grid_sizes
    % Call the function to get the numerical and exact solutions
    count = count+1;
    [U,U_exact] = FEM(n);
    errors(count) = sqrt(sum(abs(U - U_exact)^2, "all"))/(n+1)^2;
    % Display the error
    fprintf("error: %f\n", errors(count));
end

% % Loop over the grid sizes to compute errors
% for n = grid_sizes
%     % Call the function to get the numerical and exact solutions
%     [U,~] = FEM(n);
%     h = 1/(n+1);
%     sol_average = zeros(n+1,n+1);
%     exact_sol = @(x,y) (x.^2 - y.^2) .* sin(20.*x.*y);
%     for i=1:n+1
%         for j = 1:n+1
%             inter = @(x,y) (U(i,j) - exact_sol(x,y)).^2;
%             sol_average(i,j) = integral2(inter, (i-1) * h, (i) * h, (j-1) * h, (j) * h);
%         end
%     end
% 
%     % Compute the error between the numerical and exact solutions
%     error = sqrt(sum(sol_average,"all"));
% 
%     % Display the error
%    fprintf("L2 error: %f\n", error);
% end

% Define the domain and exact solution for plotting
N = 80;
x = linspace(0, 1, N+1); y = linspace(0, 1, N+1);
[X, Y] = meshgrid(x, y);
[U_num,U_exact] = FEM(N);

% Get the numerical solution and plot it
figure;
surf(X, Y, U_num);
xlabel('x');
ylabel('y');
zlabel('u_{FVM}(x, y)');
zlim([-1 1]);
title('Plot of the numerical solution using finite element method');

% Plot the original function
figure;
surf(X, Y, U_exact);
zlim([-1 1])
xlabel('x');
ylabel('y');
zlabel('u(x,y)');
title('Plot of the original function u(x,y) = (x^2 − y^2)sin(20xy)');