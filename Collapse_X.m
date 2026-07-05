function CM = Collapse_X(A)

CM = (A(:, 1:(end-1)) + A(:, 2:(end))) / 2;

