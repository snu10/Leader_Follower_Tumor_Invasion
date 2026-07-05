function [G_bd_x, G_bd_y] = BoundaryGrad(A, sc, BD_X, BD_Y)


G_bd_x = BD_X;
G_bd_y = BD_Y;

qx = diff(A, 1,2);
qy = diff(A, 1,1);

G_bd_x(:, 2:(end-1)) =qx/sc;
G_bd_y(2:(end-1), :) =qy/sc;

end