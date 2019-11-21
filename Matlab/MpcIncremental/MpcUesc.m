%Esta é uma classe para instanciar controladores MPC lineares que serão
%utilizados no minicurso de mpc da UESC.
%**************************************************************************
%                                 VARIÁVEIS
%**************************************************************************
%Propriedades
%modelo -   É um modelo em espaço de estados que define o modelo de
%           predição utilizado pelo controlador.
%**************************************************************************
%MÉTODOS
%esteMPC    -  Construto da classe
classdef MpcUesc < ModeloExpandido
   
   properties (Access = private)
      Np
      Nc
      Rbarra
      Qbarra
      F
      Phi
   end
   methods (Access = public)
      function esteMPC = MpcUesc(modelo,Np,Nc,Qbarra,Rbarra,Ts)
         NARGUMENTO = 6;
         switch nargin
            case NARGUMENTO
                argumento = {modelo Ts};
            case (NARGUMENTO-1)
                argumento = {modelo};
            otherwise
               msg = ['Erro no construtor do MPC. A quantidade de argumento'...
                  'de entrada pode ser a causa do problema.'];
               error(msg)
         end
         esteMPC@ModeloExpandido(argumento{:});
         
         
         esteMPC.Np = Np;
         if Nc>Np
            Nc = Np;
            msgbox('Não é possível atribuir Nc > Np. Nc foi alterado para Nc=Np')
         end
         esteMPC.Nc = Nc;
         esteMPC.Rbarra = Rbarra;
         esteMPC.Qbarra = Qbarra;
         esteMPC.criarMatrizesMPC()
      end % construtor
      
      %Pegar os dados da sintonia do MPC
      function sintonia = getSintonia(esteMPC)
         sintonia.horizontePredicao = esteMPC.Np;
         sintonia.horizonteControle = esteMPC.Nc;
         sintonia.MatrizRbarra = esteMPC.Rbarra;
         sintonia.MatrizQbarra = esteMPC.Qbarra;
      end
      
      function matrizes = getMatrizMpc(esteMPC)
         matrizes.F = esteMPC.F;
         matrizes.Phi = esteMPC.Phi;
      end
      
      function criarMatrizesMPC(esteMPC)       
         Npp = esteMPC.Np;
         Ncc = esteMPC.Nc;
         m = esteMPC.getModeloExpandido;
         Ae = m.A;
         Be = m.B;
         Ce = m.C;
         
         dimE = esteMPC.getDimesaoExpandida();
         nx_e = dimE.nx;
         nu_e = dimE.nu;
         ny_e = dimE.ny;
         matrizF = zeros(Npp*ny_e,nx_e);
         
         aux = zeros(Npp*ny_e,nu_e);
         for i=1:Npp
            if i ==1
               matrizF(i:i*ny_e,1:nx_e) = Ce*Ae;
               aux(i:i*ny_e,1:nu_e) = Ce*Ae*Be;
            else
               matrizF(i:i*ny_e,1:nx_e) = matrizF((i-1):(i-1)*ny_e,1:nx_e)*Ae;
               aux(i:i*ny_e,1:nu_e) = matrizF((i):(i)*ny_e,1:nx_e)*Be;
            end
         end
         esteMPC.F = matrizF;
         
         aux = [Ce*Be;aux(1:(Npp-1)*ny_e,:)];
         matrizPhi = zeros(Npp,Ncc*ny_e);
         matrizPhi(:,1:nu_e) = aux;
         for i=2:(Ncc*nu_e)
            matrizPhi((ny_e+1):end,((i-1)*nu_e+1):i*nu_e) = matrizPhi((1):(end-ny_e),((i-2)*nu_e+1):(i-1)*nu_e);
         end
         esteMPC.Phi = matrizPhi;
         %
      end
      
      function saida = calcularJk(esteMPC,Uk,Yref,xk)
         m = esteMPC.getMatrizMpc();
         s = esteMPC.getSintonia();
         phi = m.Phi;
         f = m.F;
         qbar = s.MatrizQbarra;
         rbar = s.MatrizRbarra;
         
          H = phi.'*qbar*phi + rbar;
         E = Yref - f*xk;
%          cf = -2*E.'*Qbarra*phi;
%          c = E.'*Qbarra*E;
      end
   end%end method
  
   
end % da classe