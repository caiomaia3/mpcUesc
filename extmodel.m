Ac = [0 1 0; 3 0 1; 0 1 0 ];
Bc= [1; 1; 3];
Cc=[0 1 0];
Dc=zeros(1,1);
Delta_t=1;
[Ad,Bd,Cd,Dd]=c2dm(Ac,Bc,Cc,Dc,Delta_t);