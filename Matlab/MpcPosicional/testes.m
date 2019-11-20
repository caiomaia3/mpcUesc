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
p=20;
m=5;
tSimulado = 100;

nx = 2;
nu = 1;
ny = 1;
q = 100;
r = 1;
ys =0;
uk_1 = 0;
xmk = [-1;0];
umax = 1;
umin = -1;
dumax = 0.5;

mpcUescPos =@(uk_1,xmk) issmpc(p,m,nu,ny,q,r,Ad,Bd,Cd,umax,umin,dumax,ys,uk_1,xmk); %ocultar código 
grafEstados = zeros(nx,tSimulado);
grafEntradas = zeros(nu,tSimulado);

grafJk = zeros(1,tSimulado);

for kk=1:tSimulado
   %Guarda valores dos gráficos.
   grafEstados(:,kk) = xmk;
   grafEntradas(:,kk) = uk_1;
   %Realiza a otimização com valores inciais.
   [duk,~,Jk] = mpcUescPos(uk_1,xmk); 
   
   %Atualiza variáveis
   uk = uk_1+duk;
   uk_1 = uk;
   xmk = 1*Ad*xmk+Bd*uk_1;
   grafJk(:,kk) = Jk;
end

figure()
y_ref = ys*ones(1,tSimulado);
plot(y_ref,'--r')
hold on
plot(grafEstados(1,:))
grid on

grafEstados(1,end)
figure()
plot(grafEstados(2,:))
grid on

figure()
plot(grafEstados(1,:),grafEstados(2,:))
grid on

figure()
plot(grafJk)
grid on

figure()
UB = umax*ones(1,tSimulado);
plot(UB,'--r')
hold on
LB = umin*ones(1,tSimulado);
plot(UB,'--r')
stairs(grafEntradas)
grid on





%%
