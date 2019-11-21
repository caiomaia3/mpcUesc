function [Ae,Be,Ce] = criarModeloIntegralU(A,B,C)
%CRIARMODELOINTEGRALU Acrescenta integrador no modelo de entrada.
%  Esta fun��o acrescenta uma a��o integrativa ao sistema. Isto ocorre pela
%  substitui��o da vari�vel de entrada u[t] por delta u[t]. Em decorr�ncia
%  desta mudan��o, um novo estado u[k-1] � adcionado ao modelo.
nx = size(A,1);
nu = size(B,2);
ny = size(C,1);

Ae = [A              ,B;
      zeros(nu,nx)   ,eye(nu)];
Be = [B;eye(nu)];
Ce = [C, zeros(ny,nu)];
end

