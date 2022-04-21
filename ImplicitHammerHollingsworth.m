function [X,Y] = ImplicitHammerHollingsworth(F, x_0, y_0, h, RightBorder)
% Неявный 2-х стадийный метод Хаммера-Холлингсуорта

    alfa = (3)^(1/2)/6;
    A = [1/4, 1/4-alfa; 1/4 + alfa, 1/4];
    B = [1/2, 1/2];
    C = [1/2-alfa, 1/2+alfa];
    
    v = 0; X = [];
    v = v + 1;
    X(v) = (x_0);
    Y = y_0;
    
    while (x_0 < RightBorder)
        x_1 = x_0 + h;
        
        f = @(k) [k(:,1) - F(x_0+C(1)*h, y_0+h *(A(1,1)*k(:,1) + A(1,2)*k(:,2))), 
                      k(:,2) - F(x_0+C(2)*h, y_0+h*(A(2,1)*k(:,1) + A(2,2)*k(:,2)))];
        k = fsolve(f, [F(x_0,y_0), F(x_0,y_0)]);
        k1 = k(:,1); k2=k(:,2);
%         f = @(k2) k2 - F(x_0+C(2)*h, y_0+h*(A(2,1)*k1 + A(2,2)*k2));
%         k2 = fsolve(f, k2);
        y_1 = y_0 + h*(B(1)*k1 + B(2)*k2); 
        
        x_0 = x_1;
        y_0 = y_1;
        v = v + 1;
        X(v) = x_0; 
        Y = [Y, y_0];
    end
end

