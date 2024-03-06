% This includes a slight adaption to the Mie_abcd code provided by Mätzler
% to account for absorbing medium. The scattering/absorbing efficiencies in
% absorbing medium depend on the surface of integration and thus have
% different definitions.

function result = mie_single(dia, lambda, npar, nmed)
% particle diameter (dia) and wavelength (lambda) should use the same unit.
% Alternatively, these two variables can be replaced by one variable of "a".
m=npar/nmed; % any refractive index of n = n' + 1i*n"
a=pi/lambda*dia;% size paramter (vaccum)
x=nmed*a;% size paramter (host medium)
nmax=round(2+real(x)+4*real(x)^(1/3));
n=(1:nmax); nu = (n+0.5); z=m.*x; m2=m.*m;Cn=2*n+1; 
sqx= sqrt(0.5*pi./x); sqz= sqrt(0.5*pi./z);
bx = besselj(nu, x).*sqx;
bz = besselj(nu, z).*sqz;
yx = bessely(nu, x).*sqx;
hx = bx+1i*yx;
b1x=[sin(x)/x, bx(1:nmax-1)];
b1z=[sin(z)/z, bz(1:nmax-1)];
y1x=[-cos(x)/x, yx(1:nmax-1)];
h1x= b1x+1i*y1x;
ax = x.*b1x-n.*bx;
az = z.*b1z-n.*bz;
ahx= x.*h1x-n.*hx;
an = (m2.*bz.*ax-bx.*az)./(m2.*bz.*ahx-hx.*az);
bn = (bz.*ax-bx.*az)./(bz.*ahx-hx.*az);
cn = (bx.*ahx-hx.*ax)./(bz.*ahx-hx.*az);
dn = m.*(bx.*ahx-hx.*ax)./(m2.*bz.*ahx-hx.*az);

% modified part begins
eta=2*imag(nmed)*a;
if eta==0
    f=1/2;
else
    f = (1+(eta-1)*exp(eta))/eta^2;% the flux intercepted by upper hemisphere
end
%en1 = Cn.*(abs(an).^2.*ahx.*conj(x.*hx)-abs(bn).^2.*conj(ahx).*(x.*hx));
%qsca = imag(sum(en1)/nmed)/real(nmed)/a^2/f;% near-field (inherent) 
en = Cn.*(abs(an).^2+abs(bn).^2);
qsca = (sum(en))/abs(nmed)^2/a^2/f*exp(-eta);% far-field (apparent), scaled 
en2 = Cn.*(abs(cn).^2.*(z.*bz).*conj(az)-abs(dn).^2.*conj(z.*bz).*az);
qabs = imag(sum(en2)/npar)/real(nmed)/a^2/f;% near-field (inherent) 
% modified part ends

c1n=n.*(n+2)./(n+1); c2n=Cn./n./(n+1);
anp=real(an); anpp=imag(an);
bnp=real(bn); bnpp=imag(bn);
n1=nmax-1;
clear g1;
g1(1:4,nmax)=[0; 0; 0; 0];
g1(1,1:n1)=anp(2:nmax);
g1(2,1:n1)=anpp(2:nmax);
g1(3,1:n1)=bnp(2:nmax);
g1(4,1:n1)=bnpp(2:nmax);
asy1=c1n.*(anp.*g1(1,:)+anpp.*g1(2,:)+bnp.*g1(3,:)+bnpp.*g1(4,:));
asy2=c2n.*(anp.*bnp+anpp.*bnpp);
%en=Cn.*(anp.*anp+anpp.*anpp+bnp.*bnp+bnpp.*bnpp);
asy=2*sum(asy1+asy2)/sum(en);

result = [qsca,qabs,asy];% which are Q_sca, Q_abs, and g.
end