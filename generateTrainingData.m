clear all
%% Set simulation operations
add_noise = 0; % training data will not have added noise
sumToOne = 1; % water and lipid concentrations must sum to 100% of simulated volume
save_on = 1; % can automate saving by setting this to 1

%% Set simulation constants
lipidMin = 0; % lipid concentration minimum
lipidMax = 1; % lipid concentration maximum
waterMin = 0; % water concentration minimum
waterMax = 1; % water concentration maximum
N = 100000; % number of simulations to be run
lamSWIR = [980 1200 1300]; % set of SWIR wavelengths
lamNIR = [900 930 970]; % set of NIR wavelengths
SDs = [7 10 13 16]; % S-D separations of interest
A_min = 0.2; % scattering amplitude minimum
A_max = 10; % scattering amplitude maximum
b_mean = 1.286; % mean for simulated distribution of scattering slope b
b_std = 0.521; % st. dev. for simulated distribution of scattering slope b
noise = 5; % if add_noise was set to 1, this would be the standard deviation of the zero-mean Gaussian noise
lambda0 = 980; % set reference wavelength for Power Law formula
saveName_SWIR = 'R_SWIR_train_MC';
saveName_NIR = 'R_NIR_train_MC';

load LUT_CW_multiDistance.mat % load Monte Carlo LUT

%% Preallocate space for variables and assemble SD separation and wavelength matrices
SDsep_rep_SWIR = reshape(repmat(SDs,length(lamSWIR),1),[length(lamSWIR)*length(SDs) 1]);
SDsep_rep_NIR = reshape(repmat(SDs,length(lamNIR),1),[length(lamNIR)*length(SDs) 1]);
lamSWIR_rep = repmat(lamSWIR,1,length(SDs))';
lamNIR_rep = repmat(lamNIR,1,length(SDs))';
vars_SWIR = [SDsep_rep_SWIR lamSWIR_rep]; % matrix of each S-D separation/wavelength combination (for SWIR)
vars_NIR = [SDsep_rep_NIR lamNIR_rep]; % matrix of each S-D separation/wavelength combination (for NIR)

water = zeros(N,1);
lipid = zeros(N,1);
A_temp = zeros(N,1);
b_temp = zeros(N,1);

%% Get extinction coefficients 
x = load('chromophores_SWIR.txt'); % load extinction spectra from text file
wv = x(:,1); % list of wavelengths
waterExt = x(:,2); % list of water extinction coefficients
lipidExt = x(:,3); % list of lipid extinction coefficients

for count=1:length(lamSWIR) % get extinction coefficients for SWIR wavelengths
    wvInd_SWIR = find(wv==lamSWIR(count));
    waterExt_SWIR(count) = waterExt(wvInd_SWIR);
    lipidExt_SWIR(count) = lipidExt(wvInd_SWIR);
end

for count=1:length(lamNIR) % get extinction coefficients for NIR wavelengths
    wvInd_NIR = find(wv==lamNIR(count));
    waterExt_NIR(count) = waterExt(wvInd_NIR);
    lipidExt_NIR(count) = lipidExt(wvInd_NIR);  
end


%% Begin simulations
for i=1:N
    disp(['Running simulation ',num2str(i),' of ',num2str(N)])
    
    mua_true(i,1:3) = 0; % initialize mua for current simulation (for while loop functionality)
    musp_true(i,1:3) = 0; % initialize musp for current simulation (for while loop functionality)
   
    while mua_true(i,1)<min(LUT.Mua(:)) || mua_true(i,2)<min(LUT.Mua(:)) || mua_true(i,3)<min(LUT.Mua(:)) % continue to re-randomize water and lipid concentration selection until resultant mua values are within bounds of LUT

    firstDraw(i,1) = randi([0 1],1,1); % Randomize whether water or lipid concentration is selected first
    
    if firstDraw(i,1) == 0 % lipid selected first
        lipid(i,1) = lipidMin + (lipidMax - lipidMin)*rand(1,1); % lipid selected from uniform distribution
        waterMax_temp = 1 - lipid(i,1); % temporary water maximum based on lipid concentration
        if waterMax_temp<waterMax
        else
            waterMax_temp = waterMax;
        end
        if sumToOne == 0 % for situation in which water and lipid do not need to sum to 100% (does not apply here)
            water(i,1) = waterMin + (waterMax_temp - waterMin)*rand(1,1);
        else % this applies here
            water(i,1) = 1 - lipid(i,1); % water concentration is simply the remaining volume fraction/percentage
        end
    else % water selected first -- same process here
        water(i,1) = waterMin + (waterMax - waterMin)*rand(1,1);
        lipidMax_temp = 1 - water(i,1);
        if lipidMax_temp<lipidMax
        else
            lipidMax_temp = lipidMax;
        end
        if sumToOne == 0
            lipid(i,1)= lipidMin + (lipidMax_temp - lipidMin)*rand(1,1);
        else
            lipid(i,1) = 1 - water(i,1);
        end
    end

    for j=1:3 % Beer's Law to compute mua at each wavelength 
        mua_true(i,j) = lipidExt_SWIR(j)*lipid(i,1) + waterExt_SWIR(j)*water(i,1); 
    end

    end % end while loop

    % Randomly select A and B scattering parameters from distribution
    A_temp(i,1) = -1; % initialize scattering amplitude (for while loop functionality)
    b_temp(i,1) = -1; % initialize scattering slope (for while loop functionality)
     
    while musp_true(i,3)<min(LUT.Musp(:)) || musp_true(i,1)>max(LUT.Musp(:)) || b_temp(i,1)<=0 || A_temp(i,1)<=0 % re-randomize scattering parameter selection until resultant musp values at all wavelengths are within LUT bounds
        A_temp(i,1) = A_min + (A_max - A_min)*rand(1,1); % randomly select scattering amplitude from uniform distribution
        b_temp(i,1) = normrnd(b_mean,b_std); % randomly select scattering slope from normal distribution
        for j=1:3 % compute musp values at each wavelength using power law formula
            musp_true(i,j) = A_temp(i,1)*(lamSWIR(j)/lambda0)^(-b_temp(i,1)); 
        end
    end
    
    constants = [water(i,1),lipid(i,1),A_temp(i,1),b_temp(i,1),lambda0]; % assemble array of simulation parameters for LUT mapping
    [R_SWIR(i,:)] = callLUT(constants,vars_SWIR,LUT); % get diffuse reflectance for all SWIR wavelength/S-D separation combinations 
    [R_NIR(i,:)] = callLUT(constants,vars_NIR,LUT);  % get diffuse reflectance for all NIR wavelength/S-D separation combinations 
    
    if add_noise == 1 % does not apply for training data
        R_SWIR(i,:) = R_SWIR(i,:) + normrnd(0,noise/100,[1 length(R_SWIR(i,:))]).*R_SWIR(i,:);
        R_NIR(i,:) = R_NIR(i,:) + normrnd(0,noise/100,[1 length(R_NIR(i,:))]).*R_NIR(i,:);
    end    

end % end simulation loop

if save_on==1
    save(saveName_SWIR,'R_SWIR','A_temp','b_temp','water','lipid')
    save(saveName_NIR,'R_NIR','A_temp','b_temp','water','lipid')
end