%%exemplo discretização massa-mola-amortecedor
clear variables
close all
clc

Ts = 0.010;%  Priodo de amostragem

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

clear modelo A B C K M fs

Am = modeloDiscreto.A;

Bm = modeloDiscreto.B;

Cm = modeloDiscreto.C;



%%
%Simulação
nx = size(Am,1);
tSimulacao = 200;

%Valor inicial dos estados



XHATINICIAL = [2;
          0];
       
XINICIAL = [1;
       0];
    
AMPRUIDO = 0.1;
    
xmkHat = XHATINICIAL;
xmk = XINICIAL;

%inicialização das variáveis de gráficos
grafXhat = zeros(nx,tSimulacao);
grafX = zeros(nx,tSimulacao);
kk = 1:tSimulacao;
u = sin(kk);


for kk = 1:tSimulacao
   grafXhat(:,kk) = xmkHat;
   grafX(:,kk) = xmk;
   
   xmk = Am*xmk + Bm*u(kk);
   xmkHat = Am*xmkHat + Bm*u(kk);
   
end

plot(grafXhat(1,:))
hold on
plot(grafX(1,:))

figure()

plot(grafXhat(1,:)-grafX(1,:))
grid on

%%
close all

nx = size(Am,1);

%Valor inicial dos estados
xmkHat = XHATINICIAL;
xmk = XINICIAL;

%inicialização das variáveis de gráficos
grafXhat = zeros(nx,tSimulacao);
grafX = zeros(nx,tSimulacao);
kk = 1:tSimulacao;
u = sin(kk);


for kk = 1:tSimulacao
   grafXhat(:,kk) = xmkHat;
   grafX(:,kk) = xmk;
   
   xmk = Am*xmk + Bm*u(kk)+AMPRUIDO *rand(2,1) ;
   xmkHat = Am*xmkHat + Bm*u(kk);
   
end

plot(grafXhat(1,:))
hold on
plot(grafX(1,:))

figure()

plot(grafXhat(1,:)-grafX(1,:))

figure()
plot(u)

%%
close all
plot(grafXhat(1,:)-grafX(1,:))
hold on

xmkHat = XHATINICIAL;
xmk = XINICIAL;

nx = size(Am,1);
ny = 1;

polos = 0.1*eig(Am);

Lob = place(Am',Cm',polos).';

ymk = Cm*xmk;    
    
%inicialização das variáveis de gráficos


grafXhat = zeros(nx,tSimulacao);
grafX = zeros(nx,tSimulacao);
grafY = zeros(ny,tSimulacao);

kk = 1:tSimulacao;
u = sin(0.1*kk);


for kk = 1:tSimulacao
   grafXhat(:,kk) = xmkHat;
   grafX(:,kk) = xmk;
   grafY(:,kk) = ymk;
   
   
   xmkHat = Am*xmkHat + Bm*u(kk) + Lob*(ymk-Cm*xmkHat);
   xmk = Am*xmk + Bm*u(kk) + AMPRUIDO*rand(2,1);
   ymk = Cm*xmk;

   
end

plot(grafXhat(1,:)-grafX(1,:))
grid on


figure()

plot(grafXhat(1,:),grafXhat(2,:))
hold on
plot(grafX(1,:),grafX(2,:))
% 
% figure()
% plot(grafXhat(1,:))
% hold on
% plot(grafX(1,:))
% 
% 
% figure()
% plot(u)
% 



   