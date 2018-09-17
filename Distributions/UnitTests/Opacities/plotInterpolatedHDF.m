clear all
clf %close all

fileName = './Executables/ProfileOutput.h5';

indexRadius = 180;
indexEnergy = 13;

r   = h5read( fileName, '/ProfileInfo/Radius' );
rho = h5read( fileName, '/ProfileInfo/Density' );
T   = h5read( fileName, '/ProfileInfo/Temperature' );
Ye  = h5read( fileName, '/ProfileInfo/Electron Fraction' );
E   = h5read( fileName, '/ProfileInfo/Energy' );

Ab_nue    =  h5read( fileName, ['/Opacities/','EmAb_Electron'] );
Ab_nuebar =  h5read( fileName, ['/Opacities/','EmAb_ElecAnti'] );
ES_nue    =  h5read( fileName, ['/Opacities/','ES_Electron'] );
ES_nuebar =  h5read( fileName, ['/Opacities/','ES_ElecAnti'] );

fig = figure(1);
set(gca,'FontSize',13);
subplot(2,2,1)
loglog(r, rho, r,T,'linewidth',1.5); hold on
xlabel('Radius (cm)');
ylabel('Density (g * cm-3) / Temperature (K)');
legend({'Density','Temperature'},'location','best');
plot([r(indexRadius) r(indexRadius)],...
    [min(min(rho),min(T)) max(max(rho),max(T))],...
    '-k','linewidth',1.0);

subplot(2,2,3)
semilogx(r,Ye,'linewidth',1.5); hold on
xlabel('Radius (cm)');
ylabel('Electron Fraction');
% legend('Electron Fraction','location','best');
plot([r(indexRadius) r(indexRadius)],...
    [min(Ye) max(Ye)],...
    '-k','linewidth',1.0);

subplot(2,2,2)
semilogy(E,Ab_nue(:,indexRadius),'-',...
       E,Ab_nuebar(:,indexRadius),'-',...
       E,ES_nue(:,indexRadius),'-',...
       E,ES_nuebar(:,indexRadius),'--',...
       'linewidth',1.5);hold on
xlabel('Energy (MeV)');
ylabel('1/\lambda (cm-1)');
title(['rho=',num2str(rho(indexRadius),'%.1e'),' g*cm-3',...
       ', T=', num2str(T(indexRadius),'%.1e'), 'K',...
       ', Ye=',num2str(Ye(indexRadius),'%.1e'),...
       '(black line)']);
legend({'Ab$\nu_e$',...
        'Ab $\bar{\nu}_e$',...
        'ES $\nu_e$',...
        'ES $\bar{\nu}_e$'},...
        'Interpreter','latex','location','best');   
ylim([1.e-6 1])
set(legend,'FontSize',13);

subplot(2,2,4)
loglog(r,Ab_nue(indexEnergy,:),'-',...
       r,Ab_nuebar(indexEnergy,:),'-',...
       r,ES_nue(indexEnergy,:),'-',...
       r,ES_nuebar(indexEnergy,:),'--',...
       'linewidth',1.5);hold on
xlabel('Radius (cm)');
ylabel('1/\lambda (cm-1)');
title(['Energy = ',num2str(E(indexEnergy),'%.1e'),'MeV']);
legend({'Ab $\nu_e$',...
        'Ab $\bar{\nu}_e$',...
        'ES $\nu_e$',...
        'ES $\bar{\nu}_e$'},...
        'Interpreter','latex','location','best'); 
set(legend,'FontSize',13);

x0=10;
y0=5;
width=900;
height=900;
set(gcf,'units','points','position',[x0,y0,width,height])
