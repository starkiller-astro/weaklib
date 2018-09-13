clear all
close all

fileName = './Executables/ProfileOutput.h5';

[r, rho, T, Ye, E, em_e] = ReadProfileOutputHDF...
    ( fileName, 'EmAb_Electron' );
%%% option of opacity name is
%%%    [ 'EmAb_Electron',  'EmAb_ElecAnti',
%%%      'ES_Electron',    'ES_ElecAnti']

figure(1)
subplot(2,2,1)
loglog(r, rho, r,T,'linewidth',1.5);
xlabel('Radius (cm)');
legend({'Density','Temperature'},'location','best');
subplot(2,2,3)
semilogx(r,Ye,'linewidth',1.5);
xlabel('Radius (cm)');
ylabel('Electron Fraction');
legend('Electron Fraction','location','best');
subplot(2,2,2)
loglog(E,em_e(:,100),'linewidth',1.5);hold on
xlabel('Energy (MeV)');
ylabel('Opacity');
legend(['Opacity at r=',num2str(r(100),'%E'),'cm'],'location','best');
subplot(2,2,4)
loglog(r,em_e(10,:),'linewidth',1.5);hold on
xlabel('Radius (cm)');
ylabel('Opacity');
legend(['Opacity at E=',num2str(E(10)),'MeV'],'location','best');

x0=10;
y0=5;
width=900;
height=700;
set(gcf,'units','points','position',[x0,y0,width,height])