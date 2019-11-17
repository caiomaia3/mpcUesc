function [Phi_Phi, Phi_F,Phi_R,A_e,B_e,C_e] = mpcgain(thisMPC)
%MPCGAIN - Esta função calcula as matrizes resultantaes do produtos
%mariciais nas contas do MPC.
%Como  existem algumas varáveis nas contas do MPC, cabe deixar estas contas
%em função destas. Deste modo, cabe operar o que for possivel e deixar
%estas variáveis explicitas. Assim, estas ficarão em função da entrada u0 e
%do estado inicial x0.
%% SINTAXE
%ENTRADAS
%thisMPC - um objeto instanciado da classe MpcUesc. Este objeto
%contém o modelo de predição.
%SAIDAS
%Phi_Phi - Phi*Phi
%Phi_F - Phi*F
%Phi_R - Phi*R
%A_e = Matriz A do modelo expandido
%B_e = Matriz B do modelo expandido
%C_e = Matriz C do modelo expandido


end

