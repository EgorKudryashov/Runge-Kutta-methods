% F = @(x,y)[ 2*x*(y(2))^(1/2.5)*(y(4));    %������� �� ������
%  5*x*exp(2.5*(y(3)+3))*(y(4));
%  2*x*(y(4));
%  -2*x*log(y(1))];
% x_0 = 0;
% y_0 = [1; 1; -3; 1];

%�������� ���. �������
% F = @(x,y)[ -0.04*y(1) + 10^4 * y(2) * y(3); 
%              0.04*y(1) - 10^4 * y(2) * y(3) - 3*10^7* y(2)*y(2);
%              3*10^7* y(2)*y(2)];
% x_0 = 0;
% y_0 = [1; 0; 0];
 
% %������� ������ ������� ������
% F = @(x,y) -50*(y - cos(x));
% x_0 = 0;
% y_0 = 0;

% F = @(x,y) [y(2);               %��������� ���-��� ����
%            20*((1-y(1)*y(1)) * y(2) - y(1))];
% 
% x_0 = 0;
% y_0 = [2; 0];

% % F = @(x,y) [77.27*(y(2) + y(1)*(1 - 8.375*10^(-6)* y(1) - y(2))); %����������
% %  (1/77.27) * (y(3) - (1 + y(1))* y(2));
% %   0.161* (y(1)-y(3))];
% % x_0 = 0;
% % y_0 = [1; 1; 1];

% f1 = @(x,y) -1/2*y(1);
% f2 = @(x,y) -50*(y(2) - cos(x));
% F={f1,f2};
% G= @(x,y) [-1/2*y(1); -50*(y(2) - cos(x))];
% x_0 = 0; y_0 = [1; 0];



RightBorder = 4;
h = 1/50;
tol = 10^(-6);

%[X,Y] = ExplicitRungeKuttaMethod (F, 4, x_0, y_0, h, RightBorder);
%[X,Y] = AutoStep_ERKMethod (F, 4, x_0, y_0, h, RightBorder, tol);
%[A,B] = DormanPrince (F, x_0, y_0, h, RightBorder, tol);
%[M,K] = ImplicitEuler (F, x_0, y_0, 0.1 , RightBorder);
%[R,V] = ImplicitMidPoint (F, x_0, y_0, 0.1 , RightBorder);

figure
hold on
grid on
plot (X,(Y(1,:)),'-o');
%plot (A,(B(1,:)),'-o');
%plot (M, (K(1,:)),'-*');
%plot (R, (V(1,:)),'--');
%plot (M, (K(2,:)),'--');
%plot (M, (K(3,:)),'--');
%plot (M, (K(4,:)),'--');
%plot (A,(B(2,:)),'-*');
%plot (X,Y(2,:),'-o');
%plot (X,Y(3,:),'o');
%plot (X,Y(4,:),'o');