function [X, Y] = ODIRK(F, x_0, y_0, h, rightBorder, tol)
    A = [1/4        0           0     0    0;
         1/2        1/4         0     0    0;
         17/50     -1/25       1/4    0    0;
         371/1360  -137/2720 15/544   1/4  0;
         25/24     -49/48    125/16   -85/12  1/4;
         ];
     B4 = [59/48  -17/96  225/32  -85/12];
     B5 = [25/24  -49/48  125/16  -85/12  1/4];
     C = [1/4  3/4  11/20  1/2  1];
     
    options = optimoptions('fsolve','Display','off');
    %Определяем начальные условия
    p = 4; %порядок сходимости
    fac = 0.9;
    v = 0; X = [];
    v = v + 1;
    X(v) = (x_0);
    Y = y_0;

    
    while (x_0 < rightBorder)
        f = @(k1) k1 - F(x_0+C(1)*h, y_0 + h*A(1,1)*k1);
        k1 = fsolve (f, F(x_0,y_0), options);
        f = @(k2) k2 - F(x_0+C(2)*h, y_0 + h*(A(2,1)*k1+A(2,2)*k2));
        k2 = fsolve (f, F(x_0,y_0), options);
        f = @(k3) k3 - F(x_0+C(3)*h, y_0 + h*(A(3,1)*k1 + A(3,2)*k2 + A(3,3)*k3));
        k3 = fsolve (f, F(x_0,y_0), options);
        f = @(k4) k4 - F(x_0+C(4)*h, y_0 + h*(A(4,1)*k1 + A(4,2)*k2 + A(4,3)*k3 + A(4,4)*k4));
        k4 = fsolve (f, F(x_0,y_0), options);
        f = @(k5) k5 - F(x_0+C(5)*h, y_0 + h*(A(5,1)*k1 + A(5,2)*k2 + A(5,3)*k3 + A(5,4)*k4) + A(5,5)*k5);
        k5 = fsolve (f, F(x_0,y_0), options);
        
        Y_4 = y_0 + h*(B4(1)*k1 + B4(2)*k2 + B4(3)*k3 + B4(4)*k4); 
        Y_5 = y_0 + h*(B5(1)*k1 + B5(2)*k2 + B5(3)*k3 + B5(4)*k4 + B5(5)*k5);
        
        x_0 = x_0 + h;
        y_0 = Y_5;
        
        %%%%Находим погрешность
        err = norm(abs(Y_5 - Y_4));
        
        %%%%Вычисляем следующий шаг%%%%%%%%%%%%%%%%%%%%%%%%%% 
        if (err ~= 0)
            h = fac * h * (tol/err)^(1/p);
        end
        
        v = v + 1;
        X(v) = (x_0);
        Y = [Y, y_0];   
    end

end

