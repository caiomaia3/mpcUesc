%Modelo expandido

%%exemplo discretiza��o massa-mola-amortecedor
clear variables
close all
clc

Ts = 0.1;%  Priodo de amostragem

%par�metros do sistema
K = 4.77;	%Constante da mola
M = 0.1;    %Massa
fs = 0.6;   %Constante de atrito

A =   [0     1;
      -K/M  -fs/M];
B = [0;
      1/M];  
C = [1 0];

modelo = ss(A,B,C,0); %Espa�o de estados cont�nuo
modeloDiscreto = c2d(modelo,Ts); %Discretiza��o do espa�o de estados

Ad = modeloDiscreto.A;
Bd = modeloDiscreto.B;
Cd = modeloDiscreto.C;

%%

obj = ModeloExpandido(modelo,Ts);

mExpandido = obj.getModeloExpandido;

eig(modeloDiscreto)
eig(mExpandido)

step(modeloDiscreto)
hold on
impulse(mExpandido)





