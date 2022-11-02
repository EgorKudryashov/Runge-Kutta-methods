function res = Chebyshev1Diff(s, point)
    syms x;
%     T_diff = (s * sin(s * acos(x)) / sqrt(-x^2+1));
%     res = limit (T_diff, x, point);
    t_0 = 1;
    t_1 = x;
    T = x + 1;
    for i=2:s
        T = 2*x*t_1 - t_0;
        t_0 = t_1;
        t_1 = T;
    end
    syms T_diff(x);
    T_diff(x) = diff(T, x);
    res = T_diff(point);
end

