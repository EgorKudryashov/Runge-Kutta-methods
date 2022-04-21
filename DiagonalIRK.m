function [X,Y] = DiagonalIRK(F, x_0, y_0, h, RightBorder)
% Диагональный 2-х стадийный неявный метод

    gamma = (3+(3)^(1/2))/6;
    A = [gamma 0; 1-2*gamma gamma];
    B = [1/2 1/2];
    C = [gamma 1-gamma];
    
    v = 0; X = [];
    v = v + 1;
    X(v) = (x_0);
    Y = y_0;
    
    while (x_0 < RightBorder)
        x_1 = x_0 + h;
        f = @(k1) k1 - F(x_0+C(1)*h, y_0+h*A(1,1)*k1);
        k1 = fsolve(f, F(x_0,y_0));
        f = @(k2) k2 - F(x_0+C(2)*h, y_0+h*(A(2,1)*k1 + A(2,2)*k2));
        k2 = fsolve(f, F(x_0,y_0));
        y_1 = y_0 + h*(B(1)*k1 + B(2)*k2); 
        
        x_0 = x_1;
        y_0 = y_1;
        v = v + 1;
        X(v) = x_0; 
        Y = [Y, y_0];
    end
end

