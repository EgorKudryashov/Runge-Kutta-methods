function res = Chebyshev1(s,x)
    t_0 = 1;
    t_1 = x;
    res = x + 1;
    if (s == 0)
        res = t_0;
    elseif (s==1)
            res = t_1;
    
    elseif (s > 1) 
        for i=2:s
            res = 2*x*t_1 - t_0;
            t_0 = t_1;
            t_1 = res;
        end
    end

end

