close all
clear

numerador = 2;
denominador = [1 -0.7];
Ts = 3;
Tref = 9;
Hp = 2;

Gz = tf(numerador,denominador,Ts);

uk_anterior = 0.3;
yk_anterior = 2;

y0 =2;
%t = [0,3,6,9];
[y,t] = step(0.1432*Gz);
y = y+y0;
stairs(t(1:4),y(1:4))
hold on
plot(t(1:3),y(1:3),'s')

y_dif = idpoly(denominador,numerador,'NoiseVariance',0)


nump = 1;
denp = [1,-1.4,0.45];
plant = tf(nump,denp,Ts);

plant = tf(plant);



%%




%%
%Simulação do sistema massa mola.

K = 4.77; %Constante da mola. Foi retirado de algum experimento de físca II
M = 0.1; %Massa em gramas
fs = 0.6;

A = [0,1;-K/M,-fs/M];
B = [0;1/M];
C = [1,0];
D = 0;

x0 = [0.1;0];
planta = ss(A,B,C,D);
planta_d = c2d(planta,0.1);
Ad = planta_d.A;
Bd = planta_d.B;
Cd = planta_d.C;
Dd = planta_d.D;

step(planta_d);

%%
%*************************************************************************
                     %ESTENDENDO O ESPAÇO DE ESTADOS
%*************************************************************************
nx = size(Ad,1);
nu = size(Bd,2);
ny = size(Cd,1);

Ae = [Ad, zeros(nx,ny);
   Cd*Ad, eye(ny)];
Be = [Bd;Cd*Bd];
Ce = [zeros(ny,nx), eye(ny)];
%Novas dimensões
nx_e = size(Ae,1);
nu_e = size(Be,2);
ny_e = size(Ce,1);

Np = 5;
Nc = 3;

F = zeros(Np*ny_e,nx_e);
aux = zeros(Np*ny_e,nu_e);
for i=1:Np
   if i ==1
      F(i:i*ny_e,1:nx_e) = Ce*Ae;
      aux(i:i*ny_e,1:nu_e) = Ce*Ae*Be;
   else
      F(i:i*ny_e,1:nx_e) = F((i-1):(i-1)*ny_e,1:nx_e)*Ae;
      aux(i:i*ny_e,1:nu_e) = F((i):(i)*ny_e,1:nx_e)*Be;
   end
end

aux = [Ce*Be;aux(1:(Np-1)*ny,:)];
Phi = zeros(Np,Nc*ny_e);

Phi(:,1:nu_e) = aux;
for i=2:(Nc*nu_e)
   Phi((ny_e+1):end,((i-1)*nu_e+1):i*nu_e) = Phi((1):(end-ny_e),((i-2)*nu_e+1):(i-1)*nu_e);
end

%aux(1:(Np-1)*ny,:)

%**************************************************************************
%                          Deletar depois
%******************************   __   ************************************
%******************************   ||   ************************************
%******************************  \  /  ************************************
%******************************   \/   ************************************
%**************************************************************************
% aux=Ce;
% F2 = []; 
% for i=1:Np
%    aux = aux*Ae;
%    F2 = [F2;aux];
% end
%  size(F2)==size(F)
%**************************************************************************
 
 
%%

