function [T, Y] = Milstein(F, G, h, t_0, Y_0, rightBorder)
   
    T = zeros(1, ceil(rightBorder - t_0) / h);
    rank = length(Y_0);
    Y = zeros(rank, ceil(rightBorder - t_0) / h);

    i = 1;
    T(i) = t_0;   Y(i) = Y_0;

    Ytemp = Y_0;

    rng('default');
    m = 1;
    dW = sqrt(h/m)*randn(1, ceil(m*(rightBorder-t_0)/h));
    
    diffg = diff(sym(G), 'y');
    diffG =@(t, y) diffg;
    
    while t_0 < rightBorder
    
        Winc = sum(dW((m*(i-1)+1):(m*i)));
        Ytemp = Ytemp + h*F(t_0,Ytemp);
        Ytemp = Ytemp + G(t_0,Ytemp)*Winc + 0.5 * G(t_0, Ytemp) * diffG(t_0, Ytemp) * (Winc^2 - h);
    
        t_0 = t_0 + h;
        i = i+1;
        T(i) = t_0;
        Y(i) = Ytemp;
    end
end

