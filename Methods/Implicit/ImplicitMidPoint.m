function [X,Y] = ImplicitMidPoint(F, x_0, y_0, h, RightBorder)
    
% Неявный метод средней точки
%     A = 1/2;
%     B = 1;
%     C = 1/2;
    
    v = 0; X = [];
    v = v + 1;
    X(v) = (x_0);
    Y = y_0;
    
    options = optimoptions('fsolve','Display','off');
    
    while (x_0 < RightBorder)
        x_1 = x_0 + h;
        
        f = @(k) k - F(x_0 + h/2,y_0+(h/2)*k);
        k = fsolve(f, F(x_0,y_0), options);
        y_1 = y_0 + h*k;
        
        x_0 = x_1;
        y_0 = y_1;
        v = v + 1;
        X(v) = x_1; 
        Y = [Y, y_0];
    end
end

