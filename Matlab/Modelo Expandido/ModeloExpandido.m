classdef ModeloExpandido < handle
   %MODELOEXPANDIDO Classe que gera um modelo expandido orientado a saida
   %para possiblitar um efeito integrador.
   %   Esta classe possibilita a expanção do modelo em espaço de estaddos
   %   para u
   %[ del_x_(k+1) ]  =  [Am 0m']  *  [ del_x_(k) ]   +   [  Bm  ] del_U_(k)
   %[   y_(k+1)   ]     [CmAm 1]     [   y_(k)   ]       [ CmBm ]
   %
   %     y_(k)       =  [ 0m 1 ][ del_x_(k) ]
   %                            [   y_(k)   ]
   
   %
   properties (Access = private)
      modelo
      modeloExpandido
      dim
      dim_e
   end
   
   methods (Access = public)
      function thisModel = ModeloExpandido(m,Ts)
         %MODELOEXPANDIDO construtor da classe. Veerifica consistência das
         %entradas e atribui propriedades.
         msg = ['O modelo não está em espaço de estados. Por favor forneça'...
            'um modelo definido pela função ss().'];
         assert(isa(m,'ss'),msg)
         isDiscrete = 0~=m.Ts;
         switch nargin
            case 2
               if isDiscrete
                  modelo = d2d(m,Ts,'zoh');
               else
                  modelo = c2d(m,Ts,'zoh');
               end
            case 1
               msg = ['O modelo é contínuo e precisa de um Ts para ser '...
                  'discretizado'];
               assert(isDiscrete,msg);
               modelo = m;
            otherwise
               errormsg = ['Modelo não foi definido. Uma possível causa é'...
                  'a quandtidade de argumentos.'];
               error(errormsg)
         end
         thisModel.modelo = modelo;
         thisModel.ExpandirSS();
      end %
      
      %Método para mostrar o medelo expandido configurado.
      function m = getModeloExpandido(thisModel)
         m = thisModel.modeloExpandido;
      end
      
      function dimensao = getDimesaoExpandida(thisModel)
         dimensao = thisModel.dim_e;
      end
   end % public methods
   
   methods (Access = private)
      function setDimensoes(thisModel)
         m = thisModel.modelo;
         d.nx = size(m.A,1);
         d.nu = size(m.B,2);
         d.ny = size(m.C,1);
         thisModel.dim = d;
      end
      
      function ExpandirSS(thisModel)
         
         modeloFoiDefinido = ~isempty(thisModel.modeloExpandido);
         if ~modeloFoiDefinido
            %Rotina de expandir modelo
            dimensaoDefinida = ~isempty(thisModel.dim);
            if ~dimensaoDefinida
               thisModel.setDimensoes()
            else
               disp('Dimensão já devinida anteriormente')
            end
            
            m = thisModel.modelo;
            nx = thisModel.dim.nx;
            ny = thisModel.dim.ny;
            
            Ae = [m.A, zeros(nx,ny);
                  m.C*m.A, eye(ny)];
            Be = [m.B;
                  m.C*m.B];
            Ce = [zeros(ny,nx), eye(ny)];
            
            thisModel.modeloExpandido = ss(Ae,Be,Ce,zeros(size(Ce,1),size(Be,2)),m.Ts);
            thisModel.dim_e.nx = size(Ae,1);
            thisModel.dim_e.nu = size(Be,2);
            thisModel.dim_e.ny = size(Ce,1);
            
         else
            msg = ['O modelo foi  definido anteriormente'];
            msgbox(msg);
         end 
   end
end
end

