function Z = interpolation(tau, t_0, h, F, X, Y)
    if tau < t_0
        Z = F(tau);   
    else 
        ind = bin_search(tau, X);
        theta = (tau - X(ind)) / h;
        Z = Y(:,ind) * (1-theta) + Y(:,ind+1)*theta;
    end
end

