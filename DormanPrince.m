function [X, Y] = DormanPrince(F, x_0, y_0, h, rightBorder, tol)
    A = [0 0 0 0 0 0 0;
         1/5 0 0 0 0 0 0;
         3/40 9/40 0 0 0 0 0;
         44/45 -56/15 32/9 0 0 0 0;
         19372/6561 -25360/2187 64448/6561 -212/729 0 0 0;
         9017/3168 -355/33 46732/5247 49/176 -5103/18656 0 0;
         35/384 0 500/1113 125/192 -2187/6784 11/84 0];
    B6 = [35/384 0 500/1113 125/192 -2187/6784 11/84 0];
    B7 = [5179/57600 0 7571/16695 393/640 -92097/339200 187/2100 1/40];
    C = [0 1/5 3/10 4/5 8/9 1 1];
    
    %Определяем начальные условия
    p = 5; %порядок сходимости
    fac = 0.9;
    v = 0; X = [];
    v = v + 1;
    X(v) = (x_0);
    Y = y_0;
    
    while (x_0 < rightBorder)
        Y_6 = Step (F, A,B6,C, h, x_0, y_0);
        Y_7 = Step (F, A,B7,C, h, x_0, y_0);
        
        x_0 = x_0 + h;
        y_0 = Y_7;
        
        %%%%Находим погрешность
        err = norm(abs(Y_7 - Y_6));
        
        %%%%Вычисляем следующий шаг%%%%%%%%%%%%%%%%%%%%%%%%%% 
        if (err ~= 0)
            h = fac * h * (tol/err)^(1/p);
        end
        
        v = v + 1;
        X(v) = (x_0);
        Y = [Y, y_0];   
    end
    
end

