clear all, clc
w=1; hb=1; M=1;
s=sqrt(2*hb/(M*w));
x_sz=1100;
a=9;
x=linspace(-a,a,x_sz);
N=49;
T=2*pi/w;
her=zeros(N+1,x_sz);
phi=her;
psi_0=(exp(-((x-2*s/4).^2)/8)).*cos(2*x*pi/s);%muevele al /2enexp o al 2encos
area=trapeci(conj(psi_0).*psi_0,a,x_sz);
A=1/sqrt(area);

psi_0=A*psi_0;
figure(1)
plot(x,conj(psi_0).*psi_0)
grid on
%xlr8=trapeci(conj(psi_0).*psi_0,a,x_sz)
%%
cn=zeros(1,N+1);
for n=0:N
    phi(n+1,:)=(1/pi)^.25*exp((-x.^2)/2)*(1/(sqrt(2^n*factorial(n)))).*hermiteH(n,(sqrt(2)*x)/s);
    pre_cn=phi(n+1,:).*psi_0;
    cn(n+1)=trapeci(pre_cn,a,x_sz);
end
%%
n_dts=150;
psi_nt=zeros(N+1,x_sz,n_dts+1);
psi_t=zeros(x_sz,n_dts+1);
for j=0:n_dts
    for n=0:N
        psi_nt(n+1,:,j+1)=phi(n+1,:)*cn(n+1)*exp((0-1i)*(.5+n)*w*(j*T/n_dts));
        psi_t(:,j+1)=psi_t(:,j+1).'+psi_nt(n+1,:,j+1);
    end

    funP_ti=conj(psi_t(:,j+1)).*psi_t(:,j+1);
    figure(2)
    clf
        plot(x,funP_ti,x,conj(psi_0).*psi_0);
        ylim([0 1.2]);%%%
        xlim([-a a]);
        grid on
        drawnow
end

%%
tiempo=linspace(0,T,n_dts+1);
funP=conj(psi_t(:,:)).*psi_t(:,:);
[TIEMPO,X]=meshgrid(tiempo,x);
figure(3)
surf(TIEMPO,X,funP,'EdgeColor','none')
xlabel('Tiempo')
ylabel('X')

theta=angle(psi_t);
figure(4)
surf(TIEMPO,X,theta,'EdgeColor','none')
xlabel('Tiempo')
ylabel('X')

x_esp=zeros(1,n_dts+1);
for q=0:n_dts
    pre_x_esp=conj(psi_t(:,q+1)).*x.'.*psi_t(:,q+1);
    x_esp(q+1)=trapeci(pre_x_esp,a,x_sz);
end
figure(5)
plot(tiempo,x_esp)
grid on
xlabel('Tiempo')
ylabel('<X>')
ylim([-a a]);


function [intdef]=trapeci(fun,aa,var_sz)
dx=2*aa/var_sz;
intdef=0;
for i=1:var_sz
    intdef=intdef+(dx*fun(i));
end
intdef=intdef-(dx/2)*fun(1)-(dx/2)*fun(var_sz);
end