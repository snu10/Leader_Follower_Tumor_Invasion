function A_div = Flux_Div(A_Lx ,   A_Ly ,     sc)


X_deri = diff(A_Lx, 1, 2)/sc;

Y_deri = diff(A_Ly, 1, 1)/sc;

A_div = X_deri + Y_deri;



end