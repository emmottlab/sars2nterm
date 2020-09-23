% Inhibitor data analysis

clear
clc

path = '/Users/ed/Documents/GitHub/sars2nterm/data/';

% Load data (PFU, RNA, Cytotox) 
iDat = struct();
iDat.Sora = table2array(readtable([path , 'BjornData/inhibitor_Sora.csv'])); % ML-7 
iDat.Bafe = table2array(readtable([path , 'BjornData/inhibitor_Bafe.csv'])); % ML-7 
iDat.ML7 = table2array(readtable([path , 'BjornData/inhibitor_ML7.csv'])); % ML-7 
iDat.ML18 = table2array(readtable([path , 'BjornData/inhibitor_ML18.csv'])); % ML-7 
iDat.Bosu = table2array(readtable([path , 'BjornData/inhibitor_Bosu.csv'])); % ML-7 
iDat.Sara = table2array(readtable([path , 'BjornData/inhibitor_Sara.csv'])); % ML-7 
iDat.Dasa = table2array(readtable([path , 'BjornData/inhibitor_Dasa.csv'])); % ML-7 

% File column order: 1 micromolar drug conc, 2-4 titres, 5-7 RNA, 8-10
% cytotox

% Note LOD at 80 PFU/mL for this assay (PFU)

% For main figure: Want ML7, MLCK-18, Sorafenib, and Bafetinib

samples = {'Sora','Bafe','ML7','ML18','Bosu','Sara','Dasa'};
sampleNames = {'Sorafenib','Bafetinib','ML-7','MLCK Inhibitor Peptide 18','Bosutinib','Saracatinib','Dasatinib'};
mdlDat = struct();
ci = struct();

% Note: while cell viability in Bosutinab-treated cells was estimated as 0
% in the experimental data, this was set as 1% viable for the purposes of
% model fitting.
iDat.Bosu(1,8:10) = 1;
%%

% Plot Viral Titre data
figure
for ii = 1:4 %numel(samples) % Fig 1, Plotting the 4 main hits (samples 1-4)

% hill equation sigmoid for fitting
sigmoid=@(beta,x)beta(1)+(beta(2)-beta(1))./(1+(x/beta(3)).^beta(4));

%Generate aproximate parameters to initialise the fit. 
minResponse=log10(min(min(iDat.(samples{ii})(:,2:4))));
maxResponse=log10(median(max(iDat.(samples{ii})(:,2:4))));
midResponse=mean([minResponse maxResponse]);
minDose=min(iDat.(samples{ii})(:,1));
maxDose=max(iDat.(samples{ii})(:,1));

% nonlinear fit of the data to the hill equation sigmoid
mdl = fitnlm(iDat.(samples{ii})(:,1),mean(log10(iDat.(samples{ii})(:,2:4)),2),sigmoid,[minResponse maxResponse midResponse 1]);

ci.(samples{ii}) = coefCI(mdl,0.05); % Obtain 95% confidence limits on the coefficients.

coeffs = table2array(mdl.Coefficients(:,1));
ic50 = coeffs(3); % 95% limits in ci(3)

% Now model the CC50

   subplot(2,2,ii)
    %plot the curve
    xpoints=logspace(log10(minDose),log10(maxDose),1000);
    semilogx(xpoints,sigmoid(coeffs,xpoints),'Color',[0 0 0],'LineWidth',2)
    hold on

    %Add IC50 calc
    text(ic50,mean([coeffs(1), coeffs(2)]),[sprintf('IC_{50}=%0.2g',ic50),'\muM \rightarrow'],'FontSize',10,'HorizontalAlignment','right');

    % Plot errorbars
    errorbar(iDat.(samples{ii})(:,1) , log10( mean( iDat.(samples{ii})(:,2:4),2)  ),std(log10(iDat.(samples{ii})(:,2:4)),0,2),'k','lineStyle','none'  );
    set(gca,'xscale','log')

    % Plot mean
    scatter(iDat.(samples{ii})(:,1) , log10( mean( iDat.(samples{ii})(:,2:4),2)  ),'filled','k'  );
    % Plot individual datapoints
  
    scatter([iDat.Bafe(:,1); iDat.Bafe(:,1); iDat.Bafe(:,1)],...
        [log10(iDat.(samples{ii})(:,2));log10(iDat.(samples{ii})(:,3));log10(iDat.(samples{ii})(:,2))]...
    ,'filled','r','MarkerFaceAlpha',0.3','jitter','on','jitterAmount',0.01);

% Set axis limitd
    xlim([0.03 11])  
    ylim([1.5 6.5])
    
    
    
    xlabel('Log_1_0 \muM ')% ML-7')
    ylabel('Log_1_0 PFU/mL')
    title(sampleNames{ii}) %('ML-7')
    xticks([0.01, 0.1 , 1 , 10])
    
    mdlDat.(samples{ii}) = mdl; % Save model parameters in new mdlDat structure

    clear minResponse midResponse maxResponse minDose maxDose
    line([0.03 11],[log10(80),log10(80)],'Color','k','LineStyle',':')
set(gca,'FontSize',14)    
end

% Save figure
print([path , '/Figures/Fig_Inhibitors.pdf'],'-dpdf');
%% Now the other 3 inhibitors which are in S. Figure.

figure
t = tiledlayout(2,2)

for ii = 5:7 %numel(samples) % Fig 1, Plotting the 4 main hits


%hill equation sigmoid
sigmoid=@(beta,x)beta(1)+(beta(2)-beta(1))./(1+(x/beta(3)).^beta(4));

%Generate aproximate parameters to initialise the fit. 
minResponse=log10(min(min(iDat.(samples{ii})(:,2:4))));
maxResponse=log10(median(max(iDat.(samples{ii})(:,2:4))));
midResponse=mean([minResponse maxResponse]);
minDose=min(iDat.(samples{ii})(:,1));
maxDose=max(iDat.(samples{ii})(:,1));

%[coeffs,r,J]=nlinfit(iDat.ML7(:,1),mean(log10(iDat.ML7(:,2:4)),2),sigmoid,[minResponse maxResponse midResponse 1]);
mdl = fitnlm(iDat.(samples{ii})(:,1),mean(log10(iDat.(samples{ii})(:,2:4)),2),sigmoid,[minResponse maxResponse midResponse 1]);

ci.(samples{ii}) = coefCI(mdl,0.05); % Obtain 95% confidence limits on the coefficients.

coeffs = table2array(mdl.Coefficients(:,1));
ic50 = coeffs(3); % 95% limits in ci(3)

% Now model the CC50

nexttile
 %  subplot(1,3,ii-4)
    %plot the curve
    xpoints=logspace(log10(minDose),log10(maxDose),1000);
    semilogx(xpoints,sigmoid(coeffs,xpoints),'Color',[0 0 0],'LineWidth',2)
    hold on

    %Add IC50 calc
    text(ic50,mean([coeffs(1), coeffs(2)]),['\leftarrow',sprintf('IC_{50}=%0.2g',ic50),'\muM'],'FontSize',10,'HorizontalAlignment','left');


    errorbar(iDat.(samples{ii})(:,1) , log10( mean( iDat.(samples{ii})(:,2:4),2)  ),std(log10(iDat.(samples{ii})(:,2:4)),0,2),'k','lineStyle','none'  );
    set(gca,'xscale','log')

    scatter(iDat.(samples{ii})(:,1) , log10( mean( iDat.(samples{ii})(:,2:4),2)  ),'filled','k'  );
    xlim([0.03 11]) % Means values slightly off axis edges.
    ylim([1.5 6.5])
  
    scatter([iDat.Bafe(:,1); iDat.Bafe(:,1); iDat.Bafe(:,1)],...
        [log10(iDat.(samples{ii})(:,2));log10(iDat.(samples{ii})(:,3));log10(iDat.(samples{ii})(:,2))]...
    ,'filled','r','MarkerFaceAlpha',0.3,'jitter','on','jitterAmount',0.01);

    
 %   xlabel('Log_1_0 \muM ')% ML-7')
 %   ylabel('Log_1_0 PFU/mL')
    title(sampleNames{ii}) %('ML-7')
    xticks([0.01, 0.1 , 1 , 10])
    
    mdlDat.(samples{ii}) = mdl; % Save model parameters in new mdlDat structure

    clear minResponse midResponse maxResponse minDose maxDose
    line([0.03 11],[log10(80),log10(80)],'Color','k','LineStyle',':')
set(gca,'FontSize',14)    
end

xlabel(t,'Log_1_0 \muM','FontSize',18)
ylabel(t,'Log_1_0 PFU/mL','FontSize',18)

% Save figure
print([path , '/Figures/SFig_Inhibitors.pdf'],'-dpdf');

%% Now calculate and plot cytotox data (for everything)

sigmoid=@(beta,x)beta(1)+(beta(2)-beta(1))./(1+(x/beta(3)).^beta(4));

figure
t = tiledlayout(4,2)
for ii = 1:7 %numel(samples{ii})
    ii
    
    %Generate aproximate parameters to initialise the fit.
    minResponse=log10(min(min(iDat.(samples{ii})(:,8:10))));    
    maxResponse =log10(max(max(iDat.(samples{ii})(:,8:10))));
    midResponse=mean([minResponse maxResponse]);
    minDose=min(iDat.(samples{ii})(:,1));
    maxDose=max(iDat.(samples{ii})(:,1));

    mdl = fitnlm(iDat.(samples{ii})(:,1),mean(log10(iDat.(samples{ii})(:,8:10)),2),sigmoid,[minResponse maxResponse midResponse 1]);
    coeffs = table2array(mdl.Coefficients(:,1));
    ic50 = coeffs(3); % 95% limits in ci(3)



    
%plot the curve

 nexttile   

    xpoints=logspace(log10(minDose),log10(maxDose),1000);
    semilogx(xpoints,10.^sigmoid(coeffs,xpoints),'Color','r','LineWidth',2)

    hold on

%Add CC50 calc
% numerically: where do I first hit Y < 50:
CC = 10.^sigmoid(coeffs,xpoints);
XX = CC( 1:find( CC < 50, 1 ) ); % XX = all values upto first <50 value (so can use numel(XX) to find xpoints index)

        if ii == 7% 
            %text(xpoints(numel(XX)) , 50,['\leftarrow',sprintf('IC_{50}=%0.2g',ic50),'\muM'],'FontSize',10,'Color','r','HorizontalAlignment','left');
            text(1, 50,'CC_{50}\leq0.04\muM','FontSize',10,'Color','r','HorizontalAlignment','right')
        elseif isempty(XX) == 0
            text(xpoints(numel(XX)), 50,[sprintf('CC_{50}=%0.2g',xpoints(numel(XX))),'\muM \rightarrow'],'FontSize',10,'Color','r','HorizontalAlignment','right');
        else
            text(1, 50,'CC_{50}\geq10\muM','FontSize',10,'Color','r','HorizontalAlignment','right');
        end
                
    

    errorbar(iDat.(samples{ii})(:,1) , mean( iDat.(samples{ii})(:,8:10),2)  ,std(iDat.ML7(:,8:10),0,2),'r','lineStyle','none'  );
    set(gca,'xscale','log')

    scatter(iDat.(samples{ii})(:,1) , mean( iDat.(samples{ii})(:,8:10),2)  ,'filled','r'  );
 
    scatter([iDat.Bafe(:,1); iDat.Bafe(:,1); iDat.Bafe(:,1)],...
        [iDat.(samples{ii})(:,8);iDat.(samples{ii})(:,9);iDat.(samples{ii})(:,10)]...
    ,'filled','k','MarkerFaceAlpha',0.3,'jitter','on','jitterAmount',0.01);
    
    xlim([0.03 11]) % Means values slightly off axis edges.
    ylim([0 140])
    title(sampleNames{ii})
    xticks([0.01, 0.1 , 1 , 10])
set(gca,'FontSize',14)   
end

xlabel(t,'Log_1_0 \muM','FontSize',18)
ylabel(t,'% Cell viability','FontSize',18)

% Save figure
print([path , '/Figures/SFig_InhibitorCytotox.pdf'],'-dpdf');
