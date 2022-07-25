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
%{
figure
hSurface = surf(kxa, kyb, E1)
  set(hSurface,'FaceColor',[0 0 1],'EdgeColor','none')
hold on 
hSurface = surf(kxa, kyb, E1+0.0001)
  set(hSurface,'FaceColor',[1 0 0],'EdgeColor','none')
%}

t = 0:0.01:2*pi;
x = 0.5+0.5*cos(t);
z = sin(t);
for i=1:1000
plot3(0.5+0.5*cos(t),x*0-sin(i/1000*pi/2), sin(t)-cos(i/1000*pi/2),'r')
plot3(0.5+0.49*cos(t),x*0-sin(i/1000*pi/2), 0.99*sin(t)-cos(i/1000*pi/2),'b')
xlim([-1.5,1.5])
zlim([-2.5,2.5])
ylim([-2.5, 2.5])
pause(0.01)
hold on 
end

for i=1001:3000
plot3(0.5+0.5*cos(t),x*0-sin(i/1000*pi/2), sin(t)-cos(i/1000*pi/2),'b')
plot3(0.5+0.49*cos(t),x*0-sin(i/1000*pi/2), 0.99*sin(t)-cos(i/1000*pi/2),'r')
xlim([-1.5,1.5])
zlim([-2.5,2.5])
ylim([-2.5, 2.5])
pause(0.01)
hold on 
end

for i=3001:4000
plot3(0.5+0.5*cos(t),x*0-sin(i/1000*pi/2), sin(t)-cos(i/1000*pi/2),'r')
plot3(0.5+0.49*cos(t),x*0-sin(i/1000*pi/2), 0.99*sin(t)-cos(i/1000*pi/2),'b')
xlim([-1.5,1.5])
zlim([-2.5,2.5])
ylim([-2.5, 2.5])
pause(0.01)
hold on 
end




