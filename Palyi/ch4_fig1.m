
close all 
clear all 
dt = 1/80;
t = 0:dt:1;
f = [8*(0:dt:1/8) 1+0*((1/8+dt):dt:3/8) 4-8*(3/8+dt:dt:1/2) 0*(1/2+dt:dt:1)];
plot(f) 
N = 20
k = -pi:2*pi/N:(pi-2*pi/N);
u = f-circshift(f, round(0.5/dt));
v = 2*circshift(f, -round(0.25/dt));
w = circshift(f, round(0.25/dt));

plot(t, u, 'LineWidth',1.2)

hold on 
plot(t, v, 'LineWidth', 1.2)
plot(t, w, 'LineWidth',1.2)
xlabel('t ')
ylabel('Value')
legend('u', 'v', 'w')
set(gca, "linewidth", 1, "fontsize", 18);
grid on
figure 
plot(t, sqrt(v.^2+w.^2-2*w.*v))

figure
E = zeros(length(f), 2*N);
edge = zeros(length(f), 2*N);

for i = 1:length(f) 
A = eye(N, N);
C = eye(N-1, N-1);
C = [C; zeros(1, N-1)];
C = [zeros(N,1) C];
D = eye(N-1, N-1);
D = [zeros(1, N-1); D];
D = [D zeros(N,1)];
B = [u(i) v(i); v(i) -u(i)];
Wu = [0 0; w(i) 0];
Wd = [0 w(i); 0 0];
H = kron(A,B)+kron(C, Wu)+kron(D, Wd)-0.5*eye(2*N);
[kets, energy] = eig(H);

E(i, :) = diag(energy);
[m, id] = min(abs(E(i,:))+max(E(i,:))*(E(i,:)>0));
m
edge(i,:) = (abs(kets(:, id))).^2;
%edge2(i,:) = (abs(kets(:, idx(2)))).^2;
prob = (edge(i,1:2:end))+(edge(i,2:2:end));
%prob2 = (edge2(i,1:2:end))+(edge2(i,2:2:end));

bar(prob)
title(['t = ', num2str(t(i))])
xlabel('N (Lattice site) ')
ylabel('Probability ')
ylim([0,1])
set(gca, "linewidth", 1, "fontsize", 18);
pause(0.25)
end 

plot(t, E, 'LineWidth',1.2)
xlabel('t ')
ylabel('E')
set(gca, "linewidth", 1, "fontsize", 18);
grid on


%{    
% run this block to check symbolically if our implementation is correct
N = 5;
syms u v w
A = eye(N, N);
C = eye(N-1, N-1);
C = [C; zeros(1, N-1)];
C = [zeros(N,1) C];
D = eye(N-1, N-1);
D = [zeros(1, N-1); D];
D = [D zeros(N,1)];
B = [u v; v -u];
Wu = [0 0; w 0];
Wd = [0 w; 0 0];
H = kron(A,B)+kron(C, Wu)+kron(D, Wd);
%}

