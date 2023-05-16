close all
clear all

save_on = 0; % set to 1 to automate saving of figures
%% Load data and compute errors
% load DNN-derived water and lipid estimates from csv files (saved from Python scripts)
SWIR_rec = csvread('outputDNN_SWIR_allScat_sumToOne_n1p435_g0p7_3000eps_noise5p.csv'); 
NIR_rec = csvread('outputDNN_NIR_allScat_sumToOne_n1p435_g0p7_3000eps_noise5p.csv');

% load ground truth water and lipid concentrations from test dataset files
SWIR_true = load('R_SWIR_test_MC_snell_allScat_sumToOne_n1p435_g0p7_noise.mat');
NIR_true = load('R_NIR_test_MC_snell_allScat_sumToOne_n1p435_g0p7_noise.mat');

% load loss function data from DNN training
SWIR_loss = csvread('outputDNN_SWIR_allScat_sumToOne_n1p435_g0p7_3000eps_loss.csv');
NIR_loss = csvread('outputDNN_NIR_allScat_sumToOne_n1p435_g0p7_3000eps_loss.csv');

% convert from volume fraction to volume percentage 
water_true = 100*SWIR_true.water; % ground truth water
lipid_true = 100*SWIR_true.lipid; % ground truth lipid
water_SWIR = 100*SWIR_rec(:,1); % SWIR DNN water estimate
water_NIR = 100*NIR_rec(:,1); % NIR DNN water estimate
lipid_SWIR = 100*SWIR_rec(:,2); % SWIR DNN lipid estimate
lipid_NIR = 100*NIR_rec(:,2); % NIR DNN lipid estimate

% compute absolute errors
waterAE_SWIR = (water_SWIR-water_true); 
lipidAE_SWIR = (lipid_SWIR-lipid_true); 
waterAE_NIR = (water_NIR-water_true);
lipidAE_NIR = (lipid_NIR-lipid_true);

% compute mean and standard deviation of water and lipid estimation errors
SWIR_mean = [mean(waterAE_SWIR) mean(lipidAE_SWIR)]';
SWIR_std = [std(waterAE_SWIR) std(lipidAE_SWIR)]';
NIR_mean = [mean(waterAE_NIR) mean(lipidAE_NIR)]';
NIR_std = [std(waterAE_NIR) std(lipidAE_NIR)]';

%% Create table of results
rowNames = {'Water Fraction (%)','Lipid Fraction (%)'};
varNames = {'ErrorMean_SWIR','ErrorMean_NIR','ErrorStd_SWIR','ErrorStd_NIR'};

T = table(round(SWIR_mean,2),round(NIR_mean,2),round(SWIR_std,2),round(NIR_std,2),...
    'RowNames',rowNames,'VariableNames',varNames)

%% Create figures
f1 = figure;
set(f1,'position',[10 10 700 560])
plot(water_true,water_NIR,'.','MarkerSize',5,'color',[0.5 0.5 0.5])
hold on
plot(water_true,water_SWIR,'.','MarkerSize',5,'color',[0 0 0])
plot([0 100],[0 100],'--','LineWidth',1.5,'color',[0.8 0.8 0.8])
ax = gca;
ax.PlotBoxAspectRatio = [1,1,1];
ax.TitleFontWeight = 'bold';
ax.XLabel.String = 'True Water (%)';
ax.XLabel.FontWeight = 'bold';
ax.XLabel.FontSize = 16;
ax.YLabel.String = 'Recovered Water (%)';
ax.YLabel.FontWeight = 'bold';
ax.YLabel.FontSize = 16;
ax.FontSize = 30;
thbtext.Color = [0.5,0.5,0.5];
legend({'NIR','SWIR','Identity'},'location','northwest','FontSize',16)
ylim([-5 125])

f2 = figure;
set(f2,'position',[10 10 700 560])
plot(lipid_true,lipid_NIR,'.','MarkerSize',5,'color',[0.5 0.5 0.5])
hold on
plot(lipid_true,lipid_SWIR,'.','MarkerSize',5,'color',[0 0 0])
plot([0 100],[0 100],'--','LineWidth',1.5,'color',[0.8 0.8 0.8])
ax = gca;
ax.PlotBoxAspectRatio = [1,1,1];
ax.TitleFontWeight = 'bold';
ax.XLabel.String = 'True Lipid (%)';
ax.XLabel.FontWeight = 'bold';
ax.XLabel.FontSize = 16;
ax.YLabel.String = 'Recovered Lipid (%)';
ax.YLabel.FontWeight = 'bold';
ax.YLabel.FontSize = 16;
ax.FontSize = 30;
thbtext.Color = [0.5,0.5,0.5];
legend({'NIR','SWIR','Identity'},'location','northwest','FontSize',16)
ylim([-25 105])

f3 = figure;
set(f3,'position',[10 10 700 560])
ybar = [round(SWIR_mean,2) round(NIR_mean,2); round(SWIR_std,2) round(NIR_std,2)];
b1 = bar(ybar);
b1(1).EdgeColor = 'black';
b1(1).FaceColor = 'black';
b1(2).EdgeColor = [0.5 0.5 0.5];
b1(2).FaceColor = [0.5 0.5 0.5];
xticks([1 2 3 4])
xticklabels({'Water \mu','Lipid \mu','Water \sigma','Lipid \sigma'})
ax = gca;
ax.PlotBoxAspectRatio = [1,1,1];
ax.TitleFontWeight = 'bold';
% ax.XLabel.String = 'True Lipid (%)';
ax.XLabel.FontWeight = 'bold';
ax.XLabel.FontSize = 16;
ax.YLabel.String = 'Recovered - True (%)';
ax.YLabel.FontWeight = 'bold';
ax.YLabel.FontSize = 16;
ax.FontSize = 30;
thbtext.Color = [0.5,0.5,0.5];
legend({'SWIR','NIR'},'location','northwest','FontSize',20)
xlim([0.5 4.5])
ylim([-0.4 5.9])

f4 = figure;
set(f4,'position',[10 10 700 560])
plot(NIR_loss,'-','LineWidth',1.5,'Color',[0.5 0.5 0.5])
hold on
plot(SWIR_loss,'-','LineWidth',1.5,'Color',[0 0 0])
ax = gca;
ax.YScale = 'log';
ax.PlotBoxAspectRatio = [1,1,1];
ax.TitleFontWeight = 'bold';
ax.XLabel.String = 'Epoch';
ax.XLabel.FontWeight = 'bold';
ax.XLabel.FontSize = 16;
ax.YLabel.String = 'Loss (MSE)';
ax.YLabel.FontWeight = 'bold';
ax.YLabel.FontSize = 16;
ax.FontSize = 30;
thbtext.Color = [0.5,0.5,0.5];
legend({'NIR','SWIR'},'location','northeast','FontSize',20)
xlim([-100 3000])
ylim([0.5e-5 6e-2])
yticks([1e-5 1e-4 1e-3 1e-2])

