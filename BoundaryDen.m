function [rho_bd_x, rho_bd_y] = BoundaryDen(A, BD_X, BD_Y)

rho_bd_x = BD_X;
rho_bd_y = BD_Y;

Mx1 = A(:, 1:(end-1));
Mx2 = A(:, 2:end);

rho_mean_X = (Mx1 + Mx2)/2;

My1 = A(1:(end-1), :);
My2 = A(2:end, :);

rho_mean_Y = (My1 + My2)/2;

rho_bd_x(:, 2:(end-1)) =rho_mean_X;
rho_bd_y(2:(end-1), :) =rho_mean_Y;

end