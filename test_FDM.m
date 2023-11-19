% Define grid sizes for convergence study
grid_sizes = (10:10:100);
errors = [];
count = 0;
% Loop over the grid sizes to compute errors
for n = grid_sizes
    count = count+1;
    % Call the function to get the numerical and exact solutions
    [U, U_exact] = FDM(n, n);
    % Compute the error between the numerical and exact solutions
    errors(count) = max(abs(U(:) - U_exact(:)));
    % Display the error
    fprintf("Maximum error: %f\n", errors(count));
end

% Calculate the convergence rates
% for n = grid_sizes
%     % Call the function to get the numerical and exact solutions
%     count = count+1;
%     [U,U_exact] = FDM(n,n);
%     errors(count) = sqrt(sum((U - U_exact).^2, "all")/(n^2));
%     
%     % Display the error
%     fprintf("error: %f\n", errors(count));
%     errors(count) = log(errors(count));
% end
% x = log((1:1:count));  
% 
% plot(x,errors, '-o'); 
% xlabel('X-axis');
% ylabel('Y-axis');
% title('Points Connected by Lines');

% Define the domain and exact solution for plotting
Nx = 40; Ny = 40;
dx = 1/(Nx-1); dy = 1/(Ny-1);
x = linspace(0, 1, Nx); y = linspace(0, 1, Ny);
[X, Y] = meshgrid(x, y);
U_exact = (X.^2 - Y.^2) .* sin(20.*X.*Y);

% Plot the original function
figure;
surf(X, Y, U_exact);
zlim([-1 1])
xlabel('x');
ylabel('y');
zlabel('u(x,y)');
title('Plot of the original function u(x,y) = (x^2 − y^2)sin(20xy)');

% Get the numerical solution and plot it
[U_num, ~] = FDM(Nx, Ny);

figure;
surf(X, Y, U_num);
xlabel('x');
ylabel('y');
zlabel('u_{FDM}(x, y)');
zlim([-1 1]);
title('Plot of the numerical solution using finite difference method');