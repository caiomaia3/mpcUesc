//exemplo discretização massa-mola-amortecedor
clear variables
close all
clc

Ts = 0.1;//  Priodo de amostragem

//parâmetros do sistema
K = 4.77;	//Constante da mola
M = 0.1;    //Massa
fs = 0.6;   //Constante de atrito

A =   [0     1;
      -K/M  -fs/M];
 B = [0;
      1/M];
C = [1 0];

modelo = syslin('c',A,B,C);
modeloDiscreto = dscr(modelo,Ts);

clear A B C D modelo K M fs

Ad = modeloDiscreto.A;
Bd = modeloDiscreto.B;
Cd = modeloDiscreto.C;

//sintonia do controlador

tSimulado = 40;

uk_1 = 0;
xmk = [-10;1];
function [duk,Jk] =  mpcUescPos(modelo,uk_1,xmk)
    p=10;
    m=3;
    A = modeloDiscreto.A;
    B = modeloDiscreto.B;
    C = modeloDiscreto.C;
    nx = 2;
    nu = 1;
    ny = 1;
    q = 10;
    r = 1;
    ys =0.1;
    umax = 0.5;
    umin = 0;
    dumax = 0.1;
    [duk,Jk]=issmpc(p,m,nu,ny,q,r,A,B,C,umax,umin,dumax,ys,uk_1,xmk);
endfunction
//
mpcUescPos(modeloDiscreto,uk_1,xmk)
//
//grafEstados = zeros(nx,tSimulado);
//grafEntradas = zeros(nu,tSimulado);
//
//grafJk = zeros(1,tSimulado);
//
//for kk=1:tSimulado
//   %Guarda valores dos gráficos.
//   grafEstados(:,kk) = xmk;
//   grafEntradas(:,kk) = uk_1;
//   %Realiza a otimização com valores inciais.
//   [duk,~,Jk] = mpcUescPos(uk_1,xmk); 
//   
//   %Atualiza variáveis
//   uk = uk_1+duk;
//   uk_1 = uk;
//   xmk = Ad*xmk+Bd*uk_1;
//   grafJk(:,kk) = Jk;
//end
//
//figure()
//y_ref = ys*ones(1,tSimulado);
//plot(y_ref,'--r')
//hold on
//plot(grafEstados(1,:))
//
//figure()
//plot(grafJk)
//
//figure()
//UB = umax*ones(1,tSimulado);
//plot(UB,'--r')
//hold on
//LB = umin*ones(1,tSimulado);
//plot(UB,'--r')
//stairs(grafEntradas)
