%This is an example on how to simulate polarized THz wave propagation. It
%relies on the Meridian Plane code contributed by J. Ramella-Roman, yet the 
%built-in Mie calculator has been replaced. This new calculator generates
%different parameters compared to our Mie_single function

load('n_wat.mat');
load('n_epi.mat');
%nr=linspace(1.6129,1.51,201);ni=linspace(1.5371,45.56,201);
ff = 0.01:0.01:2;
dia=150;
samp=2:2:200;% selected frequency for simulation
dop1=zeros(1,100);dop2=dop1;dop3=dop1;
RD1=dop1;RD2=dop1;RD3=dop1;
CO1=dop1;CO2=dop1;CO3=dop1;

J=1;
tic;
for jj=samp
lambda = 3e14./(ff(jj)*1e12);

% The imaginary part of refractive index needs to be positive
npar=conj(n_wat(jj));
nmed=conj(n_epi(jj));
% water in epidermis; the input needs manual changing

m=npar/(nmed);
a=pi/lambda*dia;
x=(nmed)*a;
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

% It's unnecessary to divide by f when dealing with radiative transfer
en1 = Cn.*(abs(an).^2+abs(bn).^2);
qsca = 2*(sum(en1))/abs(nmed)^2/a^2;% far-field, apparent
if imag(npar)==0
qabs = 0;
else
en2 = Cn.*(abs(cn).^2.*(z.*bz).*conj(az)-abs(dn).^2.*conj(z.*bz).*az);
qabs = 2*imag(sum(en2)/npar)/real(nmed)/a^2;% near-field, inherent
end

c1n=n.*(n+2)./(n+1); c2n=Cn./n./(n+1);
% g is not needed for this choice of phase function
% anp=real(an); anpp=imag(an);
% bnp=real(bn); bnpp=imag(bn);
% n1=nmax-1;
% g1(1:4,nmax)=[0; 0; 0; 0]; % displaced numbers used for
% g1(1,1:n1)=anp(2:nmax); % asymmetry parameter, p. 120
% g1(2,1:n1)=anpp(2:nmax);
% g1(3,1:n1)=bnp(2:nmax);
% g1(4,1:n1)=bnpp(2:nmax);
% asy1=c1n.*(anp.*g1(1,:)+anpp.*g1(2,:)+bnp.*g1(3,:)+bnpp.*g1(4,:));
% asy2=c2n.*(anp.*bnp+anpp.*bnpp);
% en=Cn.*(anp.*anp+anpp.*anpp+bnp.*bnp+bnpp.*bnpp);
% asy=2*sum(asy1+asy2)/sum(en);

% This part is quoted from the "Sphere scattering" codes shared by G. Kevin
% Zhu on File Exchange
nangles=1001; % The steps of scattering angle should be the same as defined for iquv.exe
S1=zeros(nangles,1);S2=zeros(nangles,1);
for kk=0:(nangles-1)
theta=kk/(nangles-1)*pi;
    % aux0 denotes the expression
    % assoc_legendre(nu,1,cos(theta))/sin(theta). Here I am using a
    % recursive relation to compute aux0, which avoids the numerical
    % difficulty when theta == 0 or PI.
    aux0    = zeros(1, nmax); 
    aux0(1) = -1; 
    aux0(2) = -3*cos(theta);  
    for k = 2:(nmax-1)
        aux0(k+1) = (2*k+1)/k*cos(theta)*aux0(k) - (k+1)/k*aux0(k-1);
    end
    
    % aux1 denotes the expression
    % sin(theta)*assoc_legendre_derivative(nu,1,cos(theta)).  Here I am
    % also using a recursive relation to compute aux1 from aux0,
    % which avoids numerical difficulty when theta == 0 or PI.
    aux1    = zeros(1, nmax);
    aux1(1) = cos(theta); 
    for k = 2:nmax
        aux1(k) = (k+1)*aux0(k-1)-k*cos(theta)*aux0(k);
    end
S1(kk+1)=sum(c2n.*(an.*aux0-bn.*aux1));
S2(kk+1)=sum(c2n.*(-an.*aux1+bn.*aux0));
end
% write complex array of S1 and S2 into .txt files, which are read by
% iquv.exe later
fid=fopen('S1.mci','w');
for kk=1:1001
fprintf(fid,'%8.8f %8.8f\n',real(S1(kk)),imag(S1(kk)));
end
fclose(fid);
fid2=fopen('S2.mci','w');
for kk=1:1001
fprintf(fid,'%8.8f %8.8f\n',real(S2(kk)),imag(S2(kk)));
end
fclose(fid2);

% For validation, this part can plot S1 and S2 against theta
% figure,t=tiledlayout(1,4,'TileSpacing','tight','Padding','tight');
% ax1=nexttile;semilogy((0:100)/100*180,(abs(S1).^2+abs(S2).^2)/sum(en));xlim([0,180]);xlabel('Angle (deg.)');
% ax2=nexttile;plot((0:100)/100*180,(abs(S1).^2-abs(S2).^2)./(abs(S1).^2+abs(S2).^2));xlim([0,180]);ylim([-1,1]);xlabel('Angle (deg.)');
% ax3=nexttile;plot((0:100)/100*180,2*real(S1.*conj(S2))./(abs(S1).^2+abs(S2).^2));xlim([0,180]);ylim([-1,1]);xlabel('Angle (deg.)');
% ax4=nexttile;plot((0:100)/100*180,2*imag(S1.*conj(S2))./(abs(S1).^2+abs(S2).^2));xlim([0,180]);ylim([-1,1]);xlabel('Angle (deg.)');
% hold(ax1,'on');semilogy(ax1,(0:100)/100*180,(abs(S1).^2+abs(S2).^2)/sum(en));
% hold(ax2,'on');plot(ax2,(0:100)/100*180,(abs(S1).^2-abs(S2).^2)./(abs(S1).^2+abs(S2).^2));
% hold(ax3,'on');plot(ax3,(0:100)/100*180,2*real(S1.*conj(S2))./(abs(S1).^2+abs(S2).^2));
% hold(ax4,'on');plot(ax4,(0:100)/100*180,2*imag(S1.*conj(S2))./(abs(S1).^2+abs(S2).^2));

rho     = 1/(200)^3*(150/dia)^3;     % number density, um^-3
A       = pi*dia^2/4;                % geometrical cross-sectional area, um^2
sigma_s = qsca*A;                    % scattering cross-section, um^2
mus     = sigma_s*rho*1e4;           % scattering coefficient, cm^-1
muap    = qabs*A*rho*1e4;            % particle absorbing coefficient, cm^-1
muah = 4*pi*imag(nmed)/(lambda*1e-4);% host absorbing coefficient, cm^-1
if (mus<0)||(muap<0)||(muah<0)
    fprintf('error!');
    break
end

command=strcat('iquv.exe',{' '},num2str(mus),{' '},num2str(muap),{' '},num2str(muah));
[~,result]=system(command{1});% run iquv.exe and store the output string in "result"
Stokes=split(result);
[I1,Q1,U1,V1,I2,Q2,U2,V2]=(Stokes{[4:7,11:14]});% get the output for polarization P and R.
RD1(J)=str2double(I1);RD2(J)=str2double(I2);% Diffuse Reflectance
CO1(J)=str2double(U1);CO2(J)=str2double(V2);% Co-polarization
dop1(J)=sqrt(str2double(Q1)^2+str2double(U1)^2+str2double(V1)^2)/str2double(I1);% Linear
dop2(J)=sqrt(str2double(Q2)^2+str2double(U2)^2+str2double(V2)^2)/str2double(I2);% Circular
 
fprintf('%d or %d = %.1f\n',J,jj,toc);% Show timing
J=J+1;
end

save CO1.mat CO1;
save CO2.mat CO2;
save RD1.mat RD1;
save RD2.mat RD2;
yyaxis left
plot(ff(samp),RD1);
xlabel('Frequency (THz)')
ylabel('Difffused Reflectance')
yyaxis right
plot(ff(samp),dop1);
hold on;
yyaxis right
plot(ff(samp),dop2);
ylabel('DOP')
xlim([0.15,2]);legend('Reflectance','LP','CP');
