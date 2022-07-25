clear all 
close all 
kxa = linspace(-pi, pi, 500);
kyb = linspace(-pi, pi, 500);
dt = 0.02

kyb = linspace(2*pi/3/sqrt(3)-dt, 2*pi/3/sqrt(3)+dt,  500);
kyb = kyb';
kxa = linspace(2*pi/3-dt, 2*pi/3+dt, 600);

t0 = 2.8;
td = 0*0.05*t0;
fk = 2*cos(sqrt(3)*kyb)+4*cos(sqrt(3)/2*kyb)*cos(3/2*kxa);
E1 = t0*sqrt(3+fk)+td*fk;
corr = -3/8*((kyb-2*pi/3/sqrt(3)).^2+(kxa-2*pi/3).^2)*t0.*sin(3*atan2((kyb-2*pi/3/sqrt(3)),(kxa-2*pi/3)));
E2 = 3*t0/2*sqrt((kyb-2*pi/3/sqrt(3)).^2+(kxa-2*pi/3).^2);
figure
surf(kxa, kyb, E1, 'LineStyle', 'None')
figure
surf(kxa, kyb, E2, 'LineStyle', 'None')
figure
surf(kxa, kyb, E1-E2-corr, 'LineStyle', 'None')
hold on

%surf(kxa, kyb,  E2, 'LineStyle', 'None')

corr = -3/8*((kyb-2*pi/3/sqrt(3)).^2+(kxa-2*pi/3).^2)*t0.*sin(3*atan2((kxa-2*pi/3),(kyb-2*pi/3/sqrt(3))));

figure 
surf(kxa, kyb, E1-E2-corr, 'LineStyle', 'None')
