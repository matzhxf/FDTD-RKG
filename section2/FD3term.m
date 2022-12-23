% 1/3 average
L=4; N=128; h=2*L/N; x=-L+h:h:L-h;
D=sparse(1:N-2,2:(N-1),ones(1,(N-2))/(2*h),(N-1),(N-1))+...
    sparse(2:(N-1),1:(N-2),-ones(1,(N-2))/(2*h),(N-1),(N-1));
Dx=kron(sparse(eye(N-1)),D); 
D=sparse(1:N-2,2:(N-1),ones(1,(N-2))/(h^2),(N-1),(N-1))+...
    sparse(2:(N-1),1:(N-2),ones(1,(N-2))/(h^2),(N-1),(N-1))+...
    sparse(eye(N-1)*(-2/h^2));
Dx2=kron(speye(N-1),D);
clear D
Dy=sparse(1:(N-1)*(N-2),N:(N-1)^2,ones(1,(N-1)*(N-2))/(2*h),(N-1)^2,(N-1)^2)+...
    sparse(N:(N-1)^2,1:(N-1)*(N-2),-ones(1,(N-1)*(N-2))/(2*h),(N-1)^2,(N-1)^2);
Dy2=sparse(1:(N-1)*(N-2),N:(N-1)^2,ones(1,(N-1)*(N-2))/h^2,(N-1)^2,(N-1)^2)+...
    sparse(N:(N-1)^2,1:(N-1)*(N-2),ones(1,(N-1)*(N-2))/h^2,(N-1)^2,(N-1)^2)+...
     sparse(1:(N-1)^2,1:(N-1)^2,-2*ones(1,(N-1)^2)/h^2,(N-1)^2,(N-1)^2);
[X,Y]=meshgrid(x); xx=kron(speye(N-1),sparse(1:N-1,1:N-1,x,N-1,N-1)); 
yy=kron(sparse(1:N-1,1:N-1,x,N-1,N-1),speye(N-1));
psi0=exp(-2*(X.^2+Y.^2)); psit0=(1+1i)*psi0; 
psi0=reshape(psi0.',[],1); psit0=reshape(psit0.',[],1); 
omg=1; nd=0;
k=1/64; T=2; Nmax=round(T/k);
psi1=psi0+k*psit0+k^2/2*((Dx2+Dy2)*psi0-(1+nd*abs(psi0).^2).*psi0+2*omg*(xx*Dy-yy*Dx)*psit0...
    +omg^2*(-xx.^2*Dy2-yy.^2*Dx2+xx*Dx+yy*Dy+2*xx*yy*Dx*Dy)*psi0);
psimax=zeros(1,Nmax+1); psimax(1)=max(abs(psi0));psimax(2)=max(abs(psi1));
A=speye((N-1)^2)*(1/k^2+1/2)+(omg^2*xx.^2/3-speye((N-1)^2)/2)*Dy2+...
    (omg^2*yy.^2/3-speye((N-1)^2)/2)*Dx2+(omg/k*yy-omg^2/3*xx)*Dx+...
    -(omg/k*xx+omg^2/3*yy)*Dy-omg^2*xx.*yy*Dx*Dy*2/3;
B=-speye((N-1)^2)*(1/k^2+1/2)-(omg^2*xx.^2/3-speye((N-1)^2)/2)*Dy2-...
    (omg^2*yy.^2/3-speye((N-1)^2)/2)*Dx2+(omg/k*yy+omg^2/3*xx)*Dx+...
    -(omg/k*xx-omg^2/3*yy)*Dy+omg^2*xx.*yy*Dx*Dy*2/3;
C=omg^2*(xx.^2*Dy2+yy.^2*Dx2-xx*Dx-yy*Dy-2*xx*yy*Dx*Dy)/3;
clear Dx Dy Dx2 Dy2 xx yy x psit0
for j=2:Nmax
    psi2=A\(B*psi0+2/k^2*psi1-C*psi1);
    psi0=psi1; psi1=psi2;
    psimax(j+1)=max(abs(psi1));
%     psi=(reshape(psi2,N-1,N-1)).'; pcolor(X,Y,abs(psi)); shading interp
%     drawnow
%     j*k
end
plot(0:k:T,psimax)
hold on
%psi1=(reshape(psi1,N-1,N-1)).';



