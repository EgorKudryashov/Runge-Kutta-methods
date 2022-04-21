% F = @(x,y)[ 2*x*(y(2))^(1/2.5)*(y(4));    %Задание из отчета
%  5*x*exp(2.5*(y(3)+3))*(y(4));
%  2*x*(y(4));
%  -2*x*log(y(1))];
% x_0 = 0;
% y_0 = [1; 1; -3; 1];

%Описание хим. реакций
% F = @(x,y)[ -0.04*y(1) + 10^4 * y(2) * y(3); 
%              0.04*y(1) - 10^4 * y(2) * y(3) - 3*10^7* y(2)*y(2);
%              3*10^7* y(2)*y(2)];
% x_0 = 0;
% y_0 = [1; 0; 0];
 
% %Простой пример жесткой задачи
% F = @(x,y) -50*(y - cos(x));
% %%J = [-50*sin(x), -50];
% x_0 = 0;
% y_0 = 0;

F = @(x,y) [y(2);               %уравнение Ван-Дер Поля
           20*((1-y(1)*y(1)) * y(2) - y(1))];
J = @(x,y)[0,  0,   1;
     0, - 40*y(1)*y(2) - 20, 20 - 20*y(1)^2];
x_0 = 0;
y_0 = [2; 0];

%%%Орегонатор
% % F = @(x,y) [77.27*(y(2) + y(1)*(1 - 8.375*10^(-6)* y(1) - y(2)));
% %  (1/77.27) * (y(3) - (1 + y(1))* y(2));
% %   0.161* (y(1)-y(3))];
% % x_0 = 0;
% % y_0 = [1; 1; 1];

% F = @(x,y) -1/2*y;
% x_0=0 ; y_0=1;

% %Система из жесткой задачи и нежесткой
% F= @(x,y) [-1/2*y(1); -50*(y(2) - cos(x))];
% x_0 = 0; y_0 = [1; 0];


RightBorder = 4;
h = 1/50;
tol = 10^(-8);

%[X,Y] = ExplicitRungeKuttaMethod (F, 4, x_0, y_0, h, RightBorder);
%[X,Y] = AutoStep_ERKMethod (F, 4, x_0, y_0, h, RightBorder, tol);
[A,B] = DormanPrince (F, x_0, y_0, h, RightBorder, tol);
%[M,K] = ImplicitEuler (F, x_0, y_0, 0.025 , RightBorder);
%[R,V] = ImplicitMidPoint (F, x_0, y_0, 0.025 , RightBorder);
%[C,D] = DiagonalIRK (F, x_0, y_0, 0.025, RightBorder);
[H,I] = ImplicitHammerHollingsworth(F, x_0, y_0, 0.025, RightBorder); 

figure
hold on
grid on
%plot (X,(Y(1,:)),'-o');
%plot (A,(B(1,:)),'-o');
plot (I(1,:), (I(2,:)),'--');
%plot (C, (D(1,:)),'-o');
%plot (H, (I(1,:)),'-o');
%plot (R, (V(1,:)),'--');
%plot (H, (I(2,:)),'--');
%plot (H, (I(3,:)),'--');
%plot (H, (I(4,:)),'--');

%plot (B(1,:),(B(2,:)),'-*');
%plot (A,(B(2,:)),'-*');
%plot (A,(B(3,:)),'-*');
%plot (A,(B(4,:)),'-*');

%plot (X,Y(2,:),'-o');
%plot (X,Y(3,:),'o');
%plot (X,Y(4,:),'o');