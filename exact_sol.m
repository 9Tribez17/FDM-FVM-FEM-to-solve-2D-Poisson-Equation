clc;
clear all;

N = 100;

x = linspace(0, 1, N); y = linspace(0, 1, N);
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