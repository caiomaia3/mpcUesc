%Modelo expandido

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

modelo = ss(A,B,C,0); %Espaço de estados contínuo
modeloDiscreto = c2d(modelo,Ts); %Discretização do espaço de estados

Ad = modeloDiscreto.A;
Bd = modeloDiscreto.B;
Cd = modeloDiscreto.C;


obj = ModeloExpandido(modelo,Ts);

mEstendidoY = obj.getModeloExpandido;
[Ae,Be,Ce] = criarModeloIntegralU(Ad,Bd,Cd);
De = zeros(size(Ce,1),size(Be,2));
mEstendidoU = ss(Ae,Be,Ce,De,Ts);
clear Ae Be Ce De

impulse(mEstendidoU,'.')
hold on
impulse(mEstendidoY,'r--')

legend('M. Estendido u[k-1]','M. Estendido y[k]')
title('Resposta ao Impulso')
