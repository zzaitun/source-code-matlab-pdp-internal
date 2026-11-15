clear;clc;

%Toolbox for RK4
[t,y] = ode45(@fung,[0 120],[1.5; 0]);
%[t,y] = ode45(@vdp1,[0 120],[8846930; 677; 31; 106290]);

%plot(t,y(:,1),'-o',t,y(:,2),'-o',t,y(:,3),'-o',t,y(:,4),'-o','LineWidth',2)
plot(t,y(:,1),'-o','LineWidth',2); hold on
plot(t,y(:,2),'-o','LineWidth',1.1); hold off
xlabel('Time (s)',FontSize=14);
xticks(0:5:120);
ylabel('\theta | \omega',FontSize=14);
yticklabels(["-45";"-30";"-15";"0";"15";"30";"45"]);
legend('\theta Position of angle','\omega Velocity of angle','location','north',FontSize=14);
title('Graph of Oscillation Motion',FontSize=14,FontWeight='bold');
grid on

function dydt = fung(t,y)
%Parameter
k=0.9;
mu=0.2;
m=1.9;

dydt = [y(2);
        -(mu/m)*y(2)-(k/m)*y(1)];
end