function [IE] = MD_Calculation(  )
N=750;
T=100;
L=10;
h=0.032;
part.pos.x=zeros(T,N);
part.pos.y=zeros(T,N);
part.pos.z=zeros(T,N);
part.vel.x=zeros(T,N);
part.vel.y=zeros(T,N);
part.vel.z=zeros(T,N);
%Initial state
part.pos.x(1,:)=0.013:0.013:9.75;
part.pos.y(1,:)=0.013:0.013:9.75;
part.pos.z(1,:)=0.013:0.013:9.75;
%Force exert on each particles per direction
Fx=zeros(T,N);
Fy=zeros(T,N);
Fz=zeros(T,N);
UZ=zeros(1,T);%Potential energy of the system.
KZ=zeros(1,T);%Kinetic energy of the system.
E=zeros(1,T);%Total energy of the system.
for t=[1:T];
U=zeros(1,N);%Potential energy of each particle at each time.
K=zeros(1,N);%Kinetic energy of each particle at each time.
%Force exerts on each direction on each particle at each time.
fx=zeros(1,N);
fy=zeros(1,N);
fz=zeros(1,N);
u=zeros(1,N);%Potential energy of each particle at each time.
for i=[1:N];
    for j=[1:N];
        if i==j
            break
        end
    xt=part.pos.x(t,i)-part.pos.x(t,j);
    if xt>L/2;
        xt=round(2*xt/L)*L-xt;
    elseif xt<-L/2;
        xt=xt+round(2*xt/L)*L;
    end
    yt=part.pos.y(t,i)-part.pos.y(t,j);
    if yt>L/2
        yt=round(2*yt/L)*L-yt;
    elseif yt<-L/2;
        yt=yt+L;
    end
    zt=part.pos.z(t,i)-part.pos.z(t,j);
    if zt>L/2
        zt=L-zt;
    elseif zt<-L/2;
        zt=zt+L;
    end
    rij2=xt*xt+yt*yt+zt*zt;
    rij2i=1/rij2;
    rij8i=rij2i*rij2i*rij2i*rij2i;
    rij14i=rij2i*rij2i*rij2i*rij2i*rij2i*rij2i*rij2i;
    rij6i=rij2i*rij2i*rij2i;
    rij12i=rij6i*rij6i;
    if rij2<6.25
    fx(j)=xt*(rij14i-0.5*rij8i);
    fy(j)=yt*(rij14i-0.5*rij8i);
    fz(j)=zt*(rij14i-0.5*rij8i);
    u(j)=4*(rij12i-rij6i);
    else break
    end
    end
    Fx(t,i)=sum(fx);
    Fy(t,i)=sum(fy);
    Fz(t,i)=sum(fz);
    part.pos.x(t+1,i)=part.pos.x(t,i)+part.vel.x(t,i)*h+Fx(t,i)*h*h/2;
    part.pos.y(t+1,i)=part.pos.y(t,i)+part.vel.y(t,i)*h+Fy(t,i)*h*h/2;
    part.pos.z(t+1,i)=part.pos.z(t,i)+part.vel.z(t,i)*h+Fz(t,i)*h*h/2;
    part.vel.x(t+1,i)=part.vel.x(t,i)+Fx(t,i)*h;
    part.vel.y(t+1,i)=part.vel.y(t,i)+Fy(t,i)*h;
    part.vel.z(t+1,i)=part.vel.z(t,i)+Fz(t,i)*h;
    K(i)=(.5)*(part.vel.x(t,i)*part.vel.x(t,i)+part.vel.y(t,i)*part.vel.y(t,i)+part.vel.z(t,i)*part.vel.z(t,i));
    U(i)=sum(u); 
end
UZ(t)=sum(U);
KZ(t)=sum(K);
E(t)=UZ(t)+KZ(t);
end
IE=sum(UZ)/(T*h);%Find the internal energy
plot([1:T],E,'r');hold all
end


