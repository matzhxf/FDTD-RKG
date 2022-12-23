%CN-type I
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
omg=1; omg1=omg;omg2=omg;
nd=0;
k=1/64; T=1; Nmax=round(T/k);

psi0=4*exp(-(1.5*X.^2+1.5*Y.^2)).*(X-0.5+1i*Y);%5*(X-.5+1i*Y).*(X+.5+1i*Y).*(X+1i*(Y-.5)).*(X+1i*(Y+.5)).*exp(-(2*X.^2+2*Y.^2));%2*exp(-(1.5*X.^2+1.3*Y.^2)).*(X+1i*Y).^3./(X.^2+Y.^2+.5).^(3/2); %(1+1i)*exp(-(4*X.^2+.5*Y.^2));%
psit0=0.0*exp(-(1.5*X.^2+1.5*Y.^2));%0*(1+1i)*psi0; %

psi0=reshape(psi0.',[],1); psit0=reshape(psit0.',[],1); 
%psit0=1*(yy*Dx*psi0-xx*Dy*psi0)+psit0;% rotating frame
psi1=psi0+k*psit0+k^2/2*((Dx2+Dy2)*psi0-(1+nd*abs(psi0).^2).*psi0+2*omg1*(xx*Dy-yy*Dx)*psit0...
    -omg2^2*(xx*Dy-yy*Dx)*(xx*Dy-yy*Dx)*psi0);
%psimax=zeros(1,Nmax+1); psimax(1)=max(abs(psi0));psimax(2)=max(abs(psi1));
A=speye((N-1)^2)*(1/k^2+1/2)-.5*speye((N-1)^2)*Dy2-.5*speye((N-1)^2)*Dx2...
    +omg1/k*(yy*Dx-xx*Dy)+omg2^2/2*(xx*Dy-yy*Dx)*(xx*Dy-yy*Dx);
B=-speye((N-1)^2)*(1/k^2+1/2)+.5*speye((N-1)^2)*Dy2+.5*speye((N-1)^2)*Dx2...
    +omg1/k*(yy*Dx-xx*Dy)-omg2^2/2*(xx*Dy-yy*Dx)*(xx*Dy-yy*Dx);
clear Dx Dy Dx2 Dy2 xx yy x
for j=2:Nmax
    psi2=A\(B*psi0+(2/k^2-nd*abs(psi1).^2).*psi1);
    psi0=psi1; psi1=psi2;
    %psimax(j+1)=max(abs(psi1));
    psi=(reshape(psi2,N-1,N-1)).'; 
    %pcolor(X,Y,angle(psi)); shading interp; colorbar
     %drawnow
%     j*k
end


hh=2*L/1024; xx=[-L:hh:L-hh]; mu=[0:N/2-1,-N/2:-1]*pi/L; 
psi00=[zeros(N-1,1), psi];  psi00=[zeros(1,N);psi00];
psi00=exp(1i*(xx+L).'*mu)*fft2(psi00)*exp(1i*mu.'*(xx+L))/N/N; [X,Y]=meshgrid(xx,xx);
pcolor(X,Y,abs(psi00)); shading interp; colorbar
figure
pcolor(X,Y,angle(psi00)); shading interp; colorbar



