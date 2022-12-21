function [X, Y] = RKC1(F, state, damp, x_0, y_0, h, rightBorder)
%chebyshevT(state,x);

%Объявляем начальные условия
    v = 0; X = [];
    v = v + 1;
    X(v) = (x_0);
    Y = y_0;
    
    T = @(s,w) Chebyshev1(s,w);
    omega_0 = 1 + damp/state^2;
    omega_1 = T(state, omega_0) / Chebyshev1Diff(state, omega_0);
    
    k=zeros(state+1);
    
    while (x_0 < rightBorder)
    k(1) = y_0;
    k(2) = y_0 + h * (omega_1 / omega_0) * F(x_0, k(1));
    for i=3:(state+1)
        k(i) = 2*(T(i-2, omega_0) / T(i-1, omega_0)) * (h*omega_1*F(x_0, k(i-1)) + omega_0 * k(i-1))- (T(i-3, omega_0) / T(i-1, omega_0))* k(i-2);
    end
    
    x_0 = x_0 + h;
    y_0 = k(state+1);
    
    v = v+1;
    X(v) = x_0;
    Y = [Y, y_0];
    
    end
    
end

