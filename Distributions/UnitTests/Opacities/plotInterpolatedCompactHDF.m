clear all
close all
fig = figure(1);

flag = [1,1,1,1]; % Ab, Iso, NES, TP

indexRadius = 1;
indexEnergy = 13;

fileNames = {'InterpolatedAbOutput','InterpolatedIsoOutput',...
    'InterpolatedNESOutput','InterpolatedPairOutput'};

for ii = 1:max(size(flag))
    if (flag(ii) == 1 )
        fileName = ['./Executables/',fileNames{ii},'.h5'];
        
        r   = h5read( fileName, '/ProfileInfo/Radius' );
        rho = h5read( fileName, '/ProfileInfo/Density' );
        T   = h5read( fileName, '/ProfileInfo/Temperature' );
        Ye  = h5read( fileName, '/ProfileInfo/Electron Fraction' );
        E   = h5read( fileName, '/ProfileInfo/Energy' );
        
        switch ii
            case 1 % Abem
                variNames = {'/Opacities/EmAb_Electron',...
                    '/Opacities/EmAb_ElecAnti'};
                marker = {'-r','--r'};
            case 2 % Iso
                variNames = {'/Opacities/Iso Neutrino',...
                    '/Opacities/Iso Antineutrino'};
                marker = {'-b','--b'};
            case 3 % NES
                variNames = {'/OpacitiesIMFP/NES_Electron',...
                    '/OpacitiesIMFP/NES_ElecAnti'};
                marker = {'-g','--g'};
            case 4 % Pair
                variNames = {'/OpacitiesIMFP/TP_Electron',...
                    '/OpacitiesIMFP/TP_ElecAnti'};
                marker = {'-k','--k'};
        end
        
        for jj = 1:max(size(variNames))
            plotdata = h5read( fileName, [variNames{jj}]);
            if (ii ==4)
                loglog(E,plotdata(:,indexRadius).*1e10,marker{jj},...
                    'linewidth',1.5);
            else
                loglog(E,plotdata(:,indexRadius),marker{jj},...
                    'linewidth',1.5);
            end
            hold on
        end
        
        xlabel('Energy (MeV)');
        ylabel('1/\lambda (cm-1)');
        
    end
end

legend({'$\nu_e$ Ab',...
        '$\bar{\nu}_e$ Ab',...
        '$\nu_e$ ~Iso',...
        '$\bar{\nu}_e$ ~Iso',...
        '$\nu_e$ NES',...
        '$\bar{\nu}_e$ NES',...
        '$\nu_e$ TP$\times 10^{10}$'...
        '$\bar{\nu}_e$ TP$\times 10^{10}$'},...
    'Interpreter','latex','location','best');

title(['\rho=',num2str(rho(indexRadius),'%.1e'),' g*cm-3',...
    ', T=', num2str(T(indexRadius),'%.1e'), 'K',...
    ', Ye=',num2str(Ye(indexRadius),'%.1e'),...
    '(black line)']);
ylim([3.e-11 3.e-5])
xlim([5 200])
xticks([5,10,20,50,100,200])
set(legend,'FontSize',13);

x0=10;
y0=5;
width=400;
height=600;
set(gcf,'units','points','position',[x0,y0,width,height])
saveas(fig,'Test.png','png');
