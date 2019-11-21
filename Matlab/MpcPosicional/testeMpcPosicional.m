%%exemplo discretização massa-mola-amortecedor
clear variables
close all
clc

Ts = 0.1;%  Priodo de amostragem

%parâmetros do sistema
K = 4.77;	%Constante da mola
M = 0.1;    %Massa
fs = 0.6;   %Constante de atrito

A =   [0     1;
      -K/M  -fs/M];
 B = [0;
      1/M];
C = [1 0];

modelo = ss(A,B,C,0);
modeloDiscreto = c2d(modelo,Ts);

%clear A B C D modelo

Ad = modeloDiscreto.A;
Bd = modeloDiscreto.B;
Cd = modeloDiscreto.C;

%sintonia do controlador
p=40;
m=10;
tSimulado = 100;

nx = 2;
nu = 1;
ny = 1;

q = 10;
r = 1;

ys =5;
uk_1 = 0;
%xmk = [-1;0];

umax = 100;
umin = -100;
dumax = 100;


uss = 0;
yss = 0;
xss = [0;0];

u0 = uss;
x0 = [-10;0];
y0 = Cd*x0;

yref = 100


%mpcUescPos =@(uk_1,xmk) issmpc(p,m,nu,ny,q,r,Ad,Bd,Cd,umax,umin,dumax,ys,uk_1,xmk);

[ur,yr,Jk]=ssmpc(p,m,nu,ny,nx,tSimulado,q,r,Ad,Bd,Cd,1*Ad,Bd,Cd,umax,umin,dumax,yref,uss,yss,xss,y0,u0,x0);


Yref = yref*ones(1,tSimulado);
plot(Yref,'--r')
hold on
plot(yr)
grid on


figure()
stairs(ur)
grid on

figure()
plot(Jk)




%%
