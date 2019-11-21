function [Ae,Be,Ce] = criarModeloIntegralU(A,B,C)
%CRIARMODELOINTEGRALU Acrescenta integrador no modelo de entrada.
%  Esta função acrescenta uma ação integrativa ao sistema. Isto ocorre pela
%  substituição da variável de entrada u[t] por delta u[t]. Em decorrência
%  desta mudanção, um novo estado u[k-1] é adcionado ao modelo.
nx = size(A,1);
nu = size(B,2);
ny = size(C,1);

Ae = [A              ,B;
      zeros(nu,nx)   ,eye(nu)];
Be = [B;eye(nu)];
Ce = [C, zeros(ny,nu)];
end

