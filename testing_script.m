% F = @(x,y)[ 2*x*(y(2))^(1/2.5)*(y(4));    %Задание из отчета
%  5*x*exp(2.5*(y(3)+3))*(y(4));
%  2*x*(y(4));
%  -2*x*log(y(1))];
% x_0 = 0;
% y_0 = [1 1 -3 1];

% F = @(x,y) [77.27*(y(2) + y(1)*(1 - 8.375*10^(-6)* y(1) - y(2))); %Орегонатор
%  (1/77.27) * (y(3) - (1 + y(1))* y(2));
%   0.161* (y(1)-y(3))];
% x_0 = 0;
% y_0 = [1; 1; 1];

F = @(x,y) -50*(y - cos(x));          %Простой пример жесткой задачи
x_0 = 0;
y_0 = 0;
% F = @(x,y) -y;
% x_0 = 0;
% y_0 = 1;

RightBorder = 1.5;
h = 1/50;
tol = 10^(-6);
% % 
% x_1 = x_0 + h;
% X0 = [x_1, y_0];
% [t,varf] = fsolve(f, X0)
% y_1 = h * varf + y_0

%[M,K] = ImplicitEuler (F, x_0, y_0, h, RightBorder);
[Me,Ke] = ExplicitRungeKuttaMethod (F, 2, x_0, y_0, 1/50, RightBorder);
[A,B] = RKC1 (F, 4, 20, x_0, y_0, h, RightBorder);
%[A,B] = ExplicitRungeKuttaMethod (F, 18, x_0, y_0, 1/10, RightBorder);
%[Aa,Bb] = RKC1 (F, 4, 0.05, x_0, y_0, 1/10, RightBorder);

figure
hold on
grid on
%plot (M,(K(1,:)));
plot (Me,(Ke(1,:)),'--');
plot (A,(B(1,:)),'-o');
%plot (Aa,(Bb(1,:)),'--');
%plot (M,(K(2,:)),'-*');
%plot (M,(K(3,:)),'-*');
%plot (M,(K(4,:)),'-*');