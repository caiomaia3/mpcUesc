
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

clear A B C 


nu=size(Bd,2); %nEntrada
ny=size(Cd,1); %nSaída
nx=size(Ad,1); %nEstado


%Sintonia do controlador
horizontePredicao = 15;
horizonteControle = 10;
r =1;
q =10;
matrizRbarra = r*eye(horizonteControle);
matrizQbarra = q*eye(horizontePredicao);

controlador = MpcUesc(modelo,horizontePredicao,horizonteControle,matrizQbarra,matrizRbarra,Ts);

dim = controlador.getDimesaoExpandida()
nu=size(Bd,2); %nEntrada
ny=size(Cd,1); %nSaída
nx=size(Ad,1); %nEstado


Uk=ones(horizonteControle*nu,1);
Yref=ones(horizontePredicao*ny,1);
xk=ones(nx,1);
controlador.calcularJk(Uk,Yref,xk)

M = controlador.getMatrizMpc;
F = M.F;
Phi = M.Phi;

m = controlador.getModeloExpandido();
A = m.A;
B = m.B;

eig(Phi'*Phi+matrizRbarra)

clear M
%%

r_ki = 1;
Rs_barra = ones(horizontePredicao,1);
Rs = Rs_barra*r_ki;

x0 = [0.1;0.2;0.3];
x_ki=x0;
deltaU = ((Phi.'*Phi+matrizRbarra)^-1)*(Phi.'*Rs-Phi.'*F*x0);

Kmpc = -[1, zeros(1,horizonteControle-1)] * ((Phi.'*Phi+matrizRbarra)^-1)*Phi.'*F;
%Kmpc = -((Phi.'*Phi+matrizRbarra)^-1)*Phi.'*F;
eig(A-B*Kmpc)

%%

dim_e = controlador.getDimesaoExpandida;

nx = dim_e.nx;
nu = dim_e.nu;
ny = dim_e.ny;

modelo = controlador.getModeloExpandido;
A = 1.0*modelo.A;
B = modelo.B;
nAmostra = 70;
estados = zeros(nx,nAmostra+1);
xAnterior = x0;
estados(:,1) = xAnterior ;
entradas = zeros(nu,nAmostra+1);

for iAmostra=2:nAmostra+1
   deltaU = inv(Phi.'*Phi+matrizRbarra) * (Phi.'*Rs-Phi.'*F*estados(:,iAmostra-1));
   entradas(:,iAmostra-1) = deltaU(1);
   estados(:,iAmostra) = A*estados(:,iAmostra-1)+B*entradas(:,iAmostra-1);
end

kk = 10:(10+nAmostra);
plot(kk,estados(1,:),'-o')
hold on
plot(kk,estados(2,:),'-o')

u = cumsum(entradas)

figure()
stairs(u)

% 
% iTempo = Ts*(0:nAmostras-1);
% y=lsim(modelo,U,iTempo);
% 
% GrafSaida = subplot(2,1,1);
% hold on
% grid on
% plot(iTempo, repmat(r_ki,nAmostras,1),'--k')
% plot(iTempo,y)
% 
% GrafEntrada = subplot(2,1,2);
% plot(iTempo,deltaU)
% 
% 
% %[7.2;-6.4;0;0];
% 
% %%
% % Np = 10; Nc = 4;
% % 
% % Ae = [modelo.A, zeros(dim.nx,dim.ny);
% %    modelo.C*modelo.A, eye(dim.ny)];
% % Be = [modelo.B;
% %       modelo.C*modelo.B];
% % Ce = [zeros(dim.ny,dim.nx), eye(dim.ny)];
% % 
% % modeloExpandido = ss(Ae,Be,Ce,zeros(1),Ts);
% % clear *e
% % 
% % %Novas dimensões
% % dimExp.nx = size(modeloExpandido.A,1);
% % dimExp.nu = size(modeloExpandido.B,2);
% % dimExp.ny = size(modeloExpandido.C,1);
% % 
