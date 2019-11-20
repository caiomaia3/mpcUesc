function [duk,dukk,Jk]=issmpc(p,m,nu,ny,q,r,A,B,C,umax,umin,dumax,ys,uk_1,xmk)
  
ysp=[];
for i=1:p;
  ysp=[ysp;ys];
end

Phi=[];
tha=[];
for k=1:p
    Phi=[Phi; C*A^(k)];
    tha=[tha; C*A^(k-1)*B];
end
a=tha;
Dm=a;
for iu=1:m-1; 
    a=[zeros(ny,nu);a(1:(p-1)*ny,:)];
    Dm=[Dm a];
end
Theta=Dm; 


aux=[];
for in=1:p;
  aux=[aux q];
end
Qbar=diag(aux);

clear aux; aux=[];

Mtil=[];
Itil=[];
auxM=zeros(nu,m*nu);
for in=1:m;
    aux=[aux r];
    auxM=[eye(nu) auxM(:,1:(m-1)*nu)];
    Mtil=[Mtil;auxM];
    Itil=[Itil;eye(nu)];
end
Rbar=diag(aux);

%Matriz H
H=Theta'*Qbar*Theta+Rbar;
H=(H+H')/2;

Dumax=dumax;
Umax=umax;
Umin=umin
for i=1:m-1;
 Umax=[Umax;umax];
 Umin=[Umin;umin];
 Dumax=[Dumax;dumax];
end


el = Phi*xmk-ysp;
ct = el'*Qbar*Theta;
c = (Phi*xmk-ysp)'*Qbar*(Phi*xmk-ysp);

% Ain*x <= Bin
Ain=[Mtil;-Mtil];
Bin=[Umax-Itil*uk_1;Itil*uk_1-Umin];
options=optimoptions('quadprog','display','off');  
dukk=quadprog(H,ct,Ain,Bin,[],[],-Dumax,Dumax,[],options);
duk=dukk(1:nu); 
Jk=dukk'*H*dukk+2*ct*dukk+c; % custo

