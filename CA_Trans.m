function [QC_x, QC_y, CA_Val] = CA_Trans(CA_Field_1, CA_Field_2, t, T1, T2, CTX_cond, BD_X, BD_Y)


w2 = (t-T1)/(T2-T1);
w1 = 1-w2;

if t<= T1
    w1 = 1;
    w2 = 0;
elseif t >= T2
    w1 = 0;
    w2 = 1;
end

CA_Val1 = CA_Field_1{3};
CA_Val2 = CA_Field_2{3};

CA_Val = w1*CA_Val1 + w2*CA_Val2;

GradX1 = CA_Field_1{1};
GradX2 = CA_Field_2{1};

GradY1 = CA_Field_1{2};
GradY2 = CA_Field_2{2};

CA_GradX = w1*GradX1 + w2*GradX2;
CA_GradY = w1*GradY1 + w2*GradY2;

switch CTX_cond
    % 1: CA-pos
    % 2: Neutral
    % 3: CA-neg

    case 1

        C = 40;
        [CA_W_x, CA_W_y] = BoundaryDen( CA_Weight(CA_Val, C) , BD_X, BD_Y); % Weight of chemotaxis (Chi)

        QC_x = CA_W_x.*CA_GradX;
        QC_y =CA_W_y.*CA_GradY;

    case 2

        QC_x = CA_GradX;
        QC_y = CA_GradY;

    case 3

        C = 20;
        [CA_W_x, CA_W_y] = BoundaryDen(2.5- CA_Weight(CA_Val, C) , BD_X, BD_Y);

        QC_x = CA_W_x.*CA_GradX;
        QC_y =CA_W_y.*CA_GradY;

end