function [A,B,C] = CreateRKCMethod (s, eta)
    A = zeros(s,s);
    B = zeros(1,s);
    T = @(s, x) Chebyshev1(s,x);
    omega0 = 1 + eta/s^2;
    omega1 = Chebyshev1(s, omega0)/Chebyshev1Diff(s, omega0);
     
    A(2, 1) = omega1/omega0;
    for i=3:s
            A(i, i - 1) = 2*T(i - 2, omega0)*omega1/T(i-1, omega0);
            for j=(i - 2):-1:1
                alpha = T(i - 2, omega0)/T(i-1, omega0);
                beta = T(i - 3, omega0)/T(i-1, omega0);
                A(i, j) = 2*alpha*omega0*A(i - 1, j) - beta*A(i - 2, j);
            end
    end
    
    B(s) = 2*T(s-1, omega0)*omega1/T(s, omega0);
    for j=(s - 1):-1:1
        B(j) = 2*T(s-1, omega0)*omega0*A(s, j)/ T(s, omega0) - T(s-2, omega0)*A(s - 1, j)/ T(s, omega0);
    end
    
    C = zeros(1,s);
    for i=1:s
        C(i)=sum(A(i,:));
    end
end
