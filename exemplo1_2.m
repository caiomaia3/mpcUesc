
%%
clc
clear variables
close all
%Definição do modelo

a = 0.8; b = 0.1;

Am = a; Bm = b; Cm = 1; Dm = 0; Ts = 1;
modelo = ss(Am,Bm,Cm,Dm,1);
%clear *m
%modelo.TimeUnit = 'seconds';
modelo.Name = 'exemplo1_2';

horizontePredicao = 10;
horizonteControle = 4;
rw =0;
matrizRbarra = rw*eye(horizonteControle);

controlador = MpcUesc(modelo,horizontePredicao,horizonteControle,matrizRbarra);

M = controlador.getMatrizMpc;
F = M.F;
Phi = M.Phi;


%%

r_ki = 1;
Rs_barra = ones(horizontePredicao,1);
Rs = Rs_barra*r_ki;

x0 = [0.1;0.2];
x_ki=x0;

%%
nx = 2
nu = 1
ny = 1
modelo = controlador.getModeloExpandido;
A = 1.01*modelo.A;
B = modelo.B;
nAmostra = 10;
estados = zeros(nx,nAmostra+1);
xAnterior = x0;
estados(:,1) = xAnterior ;
entradas = zeros(nu,nAmostra+1);



for iAmostra=2:nAmostra+1
   deltaU = ((Phi.'*Phi+matrizRbarra)^-1)*(Phi.'*Rs-Phi.'*F*estados(:,iAmostra-1));
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
