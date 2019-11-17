function [Phi_Phi, Phi_F,Phi_R,A_e,B_e,C_e] = mpcgain(thisMPC)
%MPCGAIN - Esta fun��o calcula as matrizes resultantaes do produtos
%mariciais nas contas do MPC.
%Como  existem algumas var�veis nas contas do MPC, cabe deixar estas contas
%em fun��o destas. Deste modo, cabe operar o que for possivel e deixar
%estas vari�veis explicitas. Assim, estas ficar�o em fun��o da entrada u0 e
%do estado inicial x0.
%% SINTAXE
%ENTRADAS
%thisMPC - um objeto instanciado da classe MpcUesc. Este objeto
%cont�m o modelo de predi��o.
%SAIDAS
%Phi_Phi - Phi*Phi
%Phi_F - Phi*F
%Phi_R - Phi*R
%A_e = Matriz A do modelo expandido
%B_e = Matriz B do modelo expandido
%C_e = Matriz C do modelo expandido


end

