% This code uses the function mie_single to generate Fig. 3 of the paper 

ff = 0.01:0.01:2;% frequency [THz]
qsca=zeros(1,200);asy=qsca;qabs=qsca;
dia = [1000,300,150,50];% particle diameter [um]
nr=linspace(1.6129,1.51,201);ni=linspace(1.5371,45.56,201);% linear fit of adiopose's property
load('npar.mat');
load('nmed.mat');

f=figure;
set(f,'Position',[100 100 550 370]);
t=tiledlayout(2,3,'TileSpacing','tight','Padding','none');
% set(t,'Unit','inches'); 
% set(t,'InnerPosition',2*[0.13 0.11 0.775 0.815]);
ax1=nexttile(1);ax2=nexttile(2);ax3=nexttile(3);ax4=nexttile(4);ax5=nexttile(5);ax6=nexttile(6);
h=[ax1,ax2,ax3,ax4,ax5,ax6];
set(h,'xdir','reverse');
set(h,'xscale','log');
set(h,'box','on');
set(h,'fontsize',9);
set(h,'LineWidth',1);
set(h,'TickDir','in');set(h,'TickLength',[0.025,0.025]);
set(h,'XMinorTick','on')
colororder([1 0 0; 0.5 1 0; 0 1 1; 0.5 0 1]);

for ii = 1:4 % air in adipose
    for jj=1:200
        lambda = 3e14./(ff(jj)*1e12);
        npar = 1;
        %nmed = real(n_epi(jj));
        nmed = nr(jj+1) + 1i*ni(jj+1)*(lambda*1e-4)/4/pi;
        u = mie_single(dia(ii), lambda, npar, nmed);
        qsca(jj) = u(1);asy(jj) = u(3);
    end
    hold(ax1,'on');plot(ax1,3e14./(ff*1e12),qsca,'linewidth',1.5);
    hold(ax4,'on');plot(ax4,3e14./(ff*1e12),asy,'linewidth',1.5);    
end

for ii = 1:4 % air in epidermis
    for jj=1:200
        npar = 1;
        nmed = conj(n_epi(jj));
        lambda = 3e14./(ff(jj)*1e12);
        u = mie_single(dia(ii), lambda, npar, nmed);
        qsca(jj) = u(1);asy(jj) = u(3);
    end
    hold(ax2,'on');plot(ax2,3e14./(ff*1e12),qsca,'linewidth',1.5);
    hold(ax5,'on');plot(ax5,3e14./(ff*1e12),asy,'linewidth',1.5);    
end

for ii = 1:4 % water in epidermis
    for jj=1:200
        npar = conj(n_wat(jj));
        nmed = conj(n_epi(jj));
        lambda = 3e14./(ff(jj)*1e12);
        u = mie_single(dia(ii), lambda, npar, nmed);
        qsca(jj) = u(1);asy(jj) = u(3);
    end
    hold(ax3,'on');plot(ax3,3e14./(ff*1e12),qsca,'linewidth',1.5);
    hold(ax6,'on');plot(ax6,3e14./(ff*1e12),asy,'linewidth',1.5);    
end

linkaxes([ax1,ax2,ax3],'xy');linkaxes([ax4,ax5,ax6],'xy');
xticks(h,[100 1000 10000]);
%xticks([ax4,ax5,ax6],[100 1000 10000]);
%set([ax1,ax2,ax3],'XAxisLocation','TOP');
%xticks([ax1,ax2,ax3],[100 1000 10000]);
%xticklabels([ax1,ax2,ax3],{'3','0.3','0.03'});
%xlabel([ax1,ax2,ax3],'Frequency (THz)');
ylabel(ax1,'Q_s_c_a');
xlabel([ax4,ax5,ax6],'Wavelength (\mum)');
ylabel(ax4,'g');
for i = 1:6
a=legend(h(i),'1000 \mum','300 \mum','150 \mum','50 \mum','location','best','FontSize',6);
a.ItemTokenSize=[10 6];a.Color='none';a.Box='off';
end
title(ax1,'air in adipose','FontSize',10)
title(ax2,'air in epidermis','FontSize',10)
title(ax3,'water in epidermis','FontSize',10)
%print('param_v2','-dtiff','-r600')



