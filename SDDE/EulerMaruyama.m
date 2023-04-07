function [T, Y] = EulerMaruyama(F, G, h, t_0, Y_0, rightBorder)

    T = zeros(1, ceil(rightBorder - t_0) / h);
    rank = length(Y_0);
    Y = zeros(rank, ceil(rightBorder - t_0) / h);

    i = 1;
    T(i) = t_0;   Y(:,i) = Y_0;

    Ytemp = Y_0;

    %rng('default');
    rng('shuffle');
    dW = sqrt(h)*randn(1, ceil(rightBorder-t_0)/h);

    while t_0 < rightBorder
        Winc = dW(i);
%         Ytemp = Ytemp + h*F(t_0,Ytemp) + G(t_0,Ytemp)*Winc;
        Ytemp = Ytemp + h*F(t_0,Ytemp);
        Ytemp = Ytemp +  G(t_0 + h,Ytemp)*Winc;
        t_0 = t_0 + h;
        i = i+1;
        T(i) = t_0;
        Y(:,i) = Ytemp;
    end

end

