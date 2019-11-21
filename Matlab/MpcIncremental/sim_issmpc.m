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
modeloDiscreto = c2d(modelo,Ts);
clear A B C modelo

[A,B,C]=immpc(modeloDiscreto.A,modeloDiscreto.B,modeloDiscreto.C);

Ap = A;
Bp = B;
Cp = C;


nu=size(B,2); %nEntrada
ny=size(C,1); %nSaída
nx=size(A,1); %nEstado

% Sintonia do MPC
p=10;% Horizonte de predição
m=3;% Horizonte de controle
tSimulacao=250;
q=[100];% Matriz de peso da saída
r=[1];% Matriz de peso das entradas

umax=[10]';
umin=[-10]';
dumax=[5]';

% Characteristics of process 
uss = [0]';
yss = [0]';

xmk=zeros(nx,1);
xmk=[-10;0;0];
xpk=xmk;
ypk=Cp*xpk;
uk_1=uss-uss;

Kf = FKalman(ny,A,C,100);


ysp=[];
for in=1:tSimulacao
    uk(:,in)=uk_1+uss;
    yk(:,in)=ypk+yss;

    if in <= 30
        ys=[0]'; % Set-point of the outputs
    else
        ys=[10]';
    end
    
    [dukk,Vk]=issmpc(p,m,nu,ny,q,r,A,B,C,umax-uss,umin-uss,dumax,ys-yss,uk_1,xmk);
    duk=dukk(1:nu); % receding horizon
    Jk(in)=Vk; % control cost
    
    %Correction of the last control input
     xmk=A*xmk+B*duk;
     ymk=C*xmk;
  if false%in==100 
      xpk=Ap*xpk+Bp*(duk+0.2*[1]'); % inserting unmeasured disturbance into plant
%       xpk=Ap*xpk+Bp*duk;
      ypk=Cp*xpk; % plant measurement 
  else
      xpk=Ap*xpk+Bp*duk;
      ypk=Cp*xpk; % plant measurement 
  end
  
  %Correction of the last measurement
  de=ypk-ymk;
  xmk=xmk+Kf*de;
  uk_1=duk+uk_1;
  ysp=[ysp ys];
end

nc=size(yk,1);
figure(1)
for j=1:nc
    subplot(nc,1,j)
    plot(yk(j,:),'k-')
    hold on
    plot(ysp(j,:),'r--')
    xlabel('tempo nT')
    in = num2str(j);
    yrot = ['y_' in];
    ylabel(yrot)
end
% legend('PV','set-point')

nc=size(uk,1);
figure(2)
for j=1:nc
subplot(nc,1,j)
    plot(uk(j,:),'k-')
    xlabel('tempo nT')
    in = num2str(j);
    yrot = ['u_' in];
    ylabel(yrot)
end

figure(3)
plot(Jk)
xlabel('tempo nT')
ylabel('Cost function')
