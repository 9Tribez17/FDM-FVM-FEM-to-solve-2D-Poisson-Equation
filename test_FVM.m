% Define grid sizes for convergence study
grid_sizes = (10:10:100);
errors = [];
count = 0;

% Loop over the grid sizes to compute errors
for n = grid_sizes
    % Call the function to get the numerical and exact solutions
    [U,~] = FVM(n,n);
    sol_average = zeros(n,n);
    dx = 1/(n-1);
    dy = 1/(n-1);
    exact_sol = @(x,y) (x.^2 - y.^2) .* sin(20.*x.*y);

    for i=1:n
        for j = 1:n
        sol_average(i,j) = integral2(exact_sol, (i-1) * dx, (i) * dx, (j-1) * dy, (j) * dy) / (dx*dy);
        end
    end

    % Compute the error between the numerical and exact solutions
    error = sum(abs(U' - sol_average),"all") * dx * dy;

    % Display the error
   fprintf("L1 error: %f\n", error);
end

% Calculate the convergence rates
for n = grid_sizes
    % Call the function to get the numerical and exact solutions
    count = count+1;
    dx = 1/(n-1);
    dy = 1/(n-1);
    [U,U_exact] = FVM(n,n);
    errors(count) = sum(abs(U - U_exact), "all")/(n*n);
    % Display the error
     fprintf( "error: %f\n",errors(count));
%     errors(count) = log(errors(count));
end
% x = log((1:1:count));  
% 
% plot(x,errors, '-o'); 
% xlabel('X-axis');
% ylabel('Y-axis');
% title('Points Connected by Lines');

% Define the domain and exact solution for plotting
Nx = 60; Ny = 60;
dx = 1/(Nx-1); dy = 1/(Ny-1);
x = linspace(0, 1, Nx); y = linspace(0, 1, Ny);
[X, Y] = meshgrid(x, y);
U_exact = (X.^2 - Y.^2) .* sin(20.*X.*Y);

% Plot the original function
% figure;
% surf(X, Y, U_exact);
% zlim([-1 1])
% xlabel('x');
% ylabel('y');
% zlabel('u(x,y)');
% title('Plot of the original function u(x,y) = (x^2 − y^2)sin(20xy)');

% Get the numerical solution and plot it
[U_num,~] = FVM(Nx, Ny);

figure;
surf(X, Y, U_num);
xlabel('x');
ylabel('y');
zlabel('u_{FVM}(x, y)');
zlim([-1 1]);
title('Plot of the numerical solution using finite volumn method');
% Plot the original function
figure;
surf(X, Y, U_exact);
zlim([-1 1])
xlabel('x');
ylabel('y');
zlabel('u(x,y)');
title('Plot of the original function u(x,y) = (x^2 − y^2)sin(20xy)');