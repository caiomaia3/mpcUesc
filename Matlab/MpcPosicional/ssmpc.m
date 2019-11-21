function [ur,yr,Jk]=ssmpc(p,m,nu,ny,nx,nsim,q,r,A,B,C,Ap,Bp,Cp,umax,umin,dumax,yspp,uss,yss,xss,y0,u0,x0)     


xpk=x0-xss ; 
xmk=xpk    ;
ypk=y0-yss ;
uk_1=u0-uss;


Psi=[];ThA=[];
for in=1:p;
    Psi=[Psi;C*A^in];
    ThA=[ThA;C*A^(in-1)*B];
end

a=ThA;
Dm=[a];
if m >= 2
    for iu=1:m-2;
        a=[zeros(ny,nu);a(1:(p-1)*ny,:)];
        Dm=[Dm a];
    end
    b=C*B;
    Ai=eye(nx);
    for in=1:p-m;
        Ai=Ai+A^in;
        b=[b;C*Ai*B];
    end
    Theta=[Dm [zeros(ny*(m-1),nu);b]];
end
if m==1
    b=C*B;
    Ai=eye(nx);
    for in=1:p-m;
        Ai=Ai+A^in;
        b=[b;C*Ai*B];
    end
    Theta=b;
end


aux=[];
for in=1:p;
  aux=[aux q];
end
Qbar=diag(aux);

clear aux; aux=[];
for in=1:m;
  aux=[aux r];
end
Rbar=diag(aux);

M=[zeros((m-1)*nu,nu) eye(nu*(m-1));zeros(nu) zeros(nu,nu*(m-1))];
Ibar=[eye(nu);zeros(nu*(m-1),nu)];
IM=eye(nu*m)-M';

%Matriz H
H=Theta'*Qbar*Theta+IM'*Rbar*IM;
H=(H+H')/2;
%
Kf = FKalman(ny,A,C,100);


Dumax=dumax;
Umax=umax-uss;
Umin=umin-uss;
for i=1:m-1;
    Umax=[Umax;umax-uss];
    Umin=[Umin;umin-uss];
    Dumax=[Dumax;dumax];
end


for in=1:nsim    
    ur(:,in)=uk_1+uss;
    yr(:,in)=ypk+yss;
    if in<= 40
        ys=yss;
    else
        ys=yspp;
    end
    ysp=[];
    for i=1:p;
        ysp=[ysp;(ys-yss)]; 
    end
    el = Psi*xmk-ysp;
    ct = el'*Qbar*Theta-uk_1'*Ibar'*Rbar*IM;
    c = (Psi*xmk-ysp)'*Qbar*(Psi*xmk-ysp)+uk_1'*Ibar'*Rbar*Ibar*uk_1;
    

    Ain=[IM;-IM];
    Bin=[Dumax+Ibar*uk_1;Dumax-Ibar*uk_1];
    options=optimoptions('quadprog','display','off');
    
    ukk=quadprog(H,ct,Ain,Bin,[],[],Umin,Umax,[],options);
    
    
    uk=ukk(1:nu);
    Jk(in)=ukk'*H*ukk+2*ct*ukk+c;
    

    xmk=A*xmk+B*uk;
    ymk=C*xmk;
    if in>=101
        %xpk=Ap*xpk+Bp*(uk+0.1*[1 .2]');
        xpk=Ap*xpk+Bp*(uk);
        ypk=Cp*xpk;
    else
        xpk=Ap*xpk+Bp*(uk);
        ypk=Cp*xpk;
    end

    de=ypk-ymk;
    xmk=xmk+Kf*de;
    uk_1=uk;
end