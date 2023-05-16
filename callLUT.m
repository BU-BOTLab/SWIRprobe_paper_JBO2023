function [R] = callLUT(constants,vars,lut)

f_water = constants(1); % water fraction
f_lipid = constants(2); % lipid fraction
a_scat = constants(3); % scattering amplitude
b_scat = constants(4); % scattering slope
lambda0 = constants(5); % reference wavelength for power law
LUT = lut; % LUT 

SD_all = vars(:,1); % full list of S-D separations
lambda_all = vars(:,2); % full list of wavelengths

x = load('chromophores_SWIR.txt'); % again, getting extinction coefficients
wv = x(:,1);
waterExt = x(:,2);
lipidExt = x(:,3);

for i=1:size(vars,1) % loop through each S-D separation/wavelength combination
   
SD = SD_all(i); % current S-D separation
lambda = lambda_all(i); % current wavelength

if SD==7 % specify which S-D separation's matrix to reference in the LUT
    R_mat = LUT.M7;
elseif SD==10
    R_mat = LUT.M10;
elseif SD==13
    R_mat = LUT.M13;
elseif SD==16
    R_mat = LUT.M16;
end

wvInd_SWIR = find(wv==lambda); % find current wavelength in extinction spectra file
E_water = waterExt(wvInd_SWIR); % get water extinction coefficient
E_lipid = lipidExt(wvInd_SWIR); % get lipid extinction coefficient

mua_temp = f_water*E_water + f_lipid*E_lipid; % Beer's Law to get temporary mua
musp_temp = a_scat*(lambda/lambda0)^(-b_scat); % power law to get temporary musp

muaInd = find(abs(mua_temp-LUT.Mua(1,:))==min(abs(mua_temp-LUT.Mua(1,:)))); % find location of LUT mua that most closely matches simulated mua
muspInd = find(abs(musp_temp-LUT.Musp(:,1))==min(abs(musp_temp-LUT.Musp(:,1)))); % find location of LUT musp that most closely matches simulated musp

R(i,1) = R_mat(muspInd,muaInd); % get diffuse reflectance from LUT

end

end

