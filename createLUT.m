%%  bookkeeping
close all
clc
clear
clear cfg
clear cfgVec
addpath(genpath('.\mcxlab'))

%% Build up model

% simulation paramters & settings
cfg.gpuid        = '11';
cfg.autopilot    = 1;
cfg.unitinmm     = 1;       % lenght of a voxel in mm
cfg.isreflect    = 0;       % Don't reflect from outside the boundary
cfg.isrefint     = 1;       % Do consider internal reflections
cfg.maxdetphoton = 1e8;     % save up to this many photons (effects memory significantly)
cfg.issaveexit   = 1;       % Save the position and direction of detected photons (can't hurt)
cfg.issaveref    = 1;       % Also save diffuse reflectance
cfg.nphoton      = 1e9;     % number of photons

% model creation
cfg.vol = uint8(ones(60,60,100)) * 1;     % creates volume and sets it all to tissue layer
cfg.vol(:,:,1)   = 0;                     % pad a layer of 0s to get diffuse reflectance

cfg.prop = [    0,  0,  1,  1;      % air optical properties [mua mus g n]
                0,  0,  0,  0];      % tissue optical properties (initializaing with zeros)
gfactor = 0.7; % set sample anisotropy factor
nind = 1.435; % set sample refractive index 
nair = 1; % set air refractive index
nglass = 1.5; % set glass refractive index (for photodiode package window)
theta_crit = asind(nair/nglass); % computation of critical angle for total internal reflection
muspRange=linspace(0.1,10,1000); % set range of musp values to simulate
muaRange=linspace(0.001,0.2,1000); % set range of mua values to simulate
 
% source and detector 
cfg.srcpos = [18, 30, 1];           % source position (x,y,z)
cfg.srcdir = [0,  0, 1];            % source direction
cfg.detpos = [25, 30, 2, 5.3/2;     % detector positions and glass window radius (x,y,z,r)
              28, 30, 2, 5.3/2;
              31, 30, 2, 5.3/2;
              34, 30, 2, 5.3/2];    
          
dist2sens = 2.5; % distance between photodiode window and photosensitive area
radius = 1.5; % radius of photosensitive portion of photodiode

% Timing information
cfg.tstart = 0;       % Time the simulation starts
cfg.tstep  = 2e-10;    % Steps to take
cfg.tend   = 5e-9;    % When to end.  The output will have [tstart:tstep:tend] slices at each of the different time points

% Initializing LUT variables for storage
Mua = zeros(length(muspRange),length(muaRange));
Musp = zeros(length(muspRange),length(muaRange));
M7_count = zeros(length(muspRange),length(muaRange));
M10_count = zeros(length(muspRange),length(muaRange));
M13_count = zeros(length(muspRange),length(muaRange));
M16_count = zeros(length(muspRange),length(muaRange));

%% Model: creates a stucture with multiple models to run (different wavelengths and updated mua & mus)
% White monte carlo -- will scale partial path lengths to account for mua
clear cfgVec

for i = 1:length(muspRange)
    disp(['Musp ',num2str(i),' of ',num2str(length(muspRange))])

    mus_temp = muspRange(i)/(1-gfactor); % compute mus from desired musp and g
    cfgVec(i) = cfg; % set current iteration's model to predefined model  
    cfgVec(i).prop(2,:) = [0, mus_temp, gfactor, nind]; % optical properties of tissue layer (zero absorption, will scale later)
    
    [flux,det] = mcxlab(cfgVec(i)); % this runs the Monte Carlo simulation using MCX
    
    theta1 = acosd(-det.v(:,3)); % compute incident angle for each photon/glass window interface
    theta2 = asind((nind/nglass)*sind(theta1)); % Snell's law to get angle of each photon after transitioning from sample to glass window
    theta2_temp = theta2; % store these angles in temporary variable
    theta2(theta2>theta_crit) = 0; % if angle after crossing from sample to glass is greater than critical angle, temporarily set angle to zero (just to avoid imaginary numbers)
    theta3 = asind((nglass/nair)*sind(theta2)); % Snell's law to get angle of each photon after transitioning from glass to air (prior to hitting photosensitive area)
    theta3(theta2_temp>theta_crit) = 89; % reject photons that previously exceeded critical angle (those temporarily set to zero) by deflecting them nearly 90 degrees away
    theta = theta3; % this is the final array of photon angles in Z-direction after passing through glass window
    alpha = atand(det.v(:,2)./det.v(:,1)); % compute angle of photon travel in XY space
    deltaX = dist2sens*tand(theta).*cosd(alpha); % compute distance traveled in X direction between glass window and photosensitive area
    deltaY = dist2sens*tand(theta).*sind(alpha); % compute distance traveled in Y direction between glass window and photosensitive area
    posXY = [(det.p(:,1) + deltaX) (det.p(:,2) + deltaY)]; % compute updated XY positions of all photons at time of reaching photosensitive area
    
    sd7_pos = posXY(det.detid==1,:); % group detected photon XY positions by S-D separation
    sd10_pos = posXY(det.detid==2,:); 
    sd13_pos = posXY(det.detid==3,:);
    sd16_pos = posXY(det.detid==4,:);
    
    sd7_dist = vecnorm(sd7_pos - [25 30],2,2); % compute XY distance between detected photons and center of photosensitive area
    sd10_dist = vecnorm(sd10_pos - [28 30],2,2);
    sd13_dist = vecnorm(sd13_pos - [31 30],2,2);
    sd16_dist = vecnorm(sd16_pos - [34 30],2,2);

    sd_grouped_7 = det.ppath(find(det.detid==1)); % group detected photon pathlengths by S-D separation
    sd_grouped_10 = det.ppath(find(det.detid==2));
    sd_grouped_13 = det.ppath(find(det.detid==3));
    sd_grouped_16 = det.ppath(find(det.detid==4));
    
    sd_trimmed_7 = sd_grouped_7(sd7_dist <= radius); % reject "detected" photons that did not hit photosensitive area after traveling through glass window 
    sd_trimmed_10 = sd_grouped_10(sd10_dist <= radius);
    sd_trimmed_13 = sd_grouped_13(sd13_dist <= radius);
    sd_trimmed_16 = sd_grouped_16(sd16_dist <= radius);
    
    for j=1:length(muaRange) % scaling photon counts for range of simulated mua values
        Mua(i,j)=muaRange(j); % assemble 2D matrix of mua values
        Musp(i,j)=muspRange(i); % assemble 2D matrix of musp values
        M7_count(i,j)=sum(exp(-muaRange(j)*sd_trimmed_7)); % Beer-Lambert Law to scale photon weight, sum all photon weights together
        M10_count(i,j)=sum(exp(-muaRange(j)*sd_trimmed_10));
        M13_count(i,j)=sum(exp(-muaRange(j)*sd_trimmed_13));
        M16_count(i,j)=sum(exp(-muaRange(j)*sd_trimmed_16));
    end
end

LUT.Mua=Mua; % final LUT matrix of mua values
LUT.Musp=Musp; % final LUT matrix of musp values
LUT.M7 = double(M7_count)/cfg.nphoton; % divide sum of photon weights by number of launched photons to get diffuse reflectance (per S-D separation)
LUT.M10 = double(M10_count)/cfg.nphoton;
LUT.M13 = double(M13_count)/cfg.nphoton;
LUT.M16 = double(M16_count)/cfg.nphoton;

% saving data
save('LUT_CW_multiDistance','LUT','-v7.3');

