function [T, meanY] = TestEquastionROCK(alpha, beta, gamma, delta, s, eta, mode, numExp)
    %%% mode = 1 - ROCK_SDDE
    %%% mode = 2 - ROCK_SDDE2
    h = 1;
    rightBorder = 50;

    F= @(t,y,z) alpha*y + beta*z(1,1);  G = @(t,y,z) sqrt(gamma)*y + sqrt(delta)*z(1,1);
    history = @(t,y) 1;
    delay = 2;

    t_0 = 0;
    Yzero = 1;

    if (mode == 1)
        for i=1:numExp
            [T, ROCK(i,:)] = ROCK_SDDE(F, G, h, t_0, Yzero, delay, history, rightBorder, s, eta);
        end
    end
    if (mode == 2)
        for i=1:numExp
            [T, ROCK(i,:)] = ROCK_SDDE2(F, G, h, t_0, Yzero, delay, history, rightBorder, s, eta);
        end
    end
    
    for i=1:length(T)
        meanY(i) = mean(ROCK(:,i).^2);
    end
end

