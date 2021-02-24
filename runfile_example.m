%Example run for ICCS

% Kamila Mustafina

path = 'data';
imagefiles_ch1 = dir(fullfile(path, '\ch1\*.tiff'));
imagefiles_ch2 = dir(fullfile(path, '\ch2\*.tiff'));
ICS_ch1_series = {};
ICS_ch2_series = {};
ICCS_series = {};
filenames = {};

%importing pictures

for ii = 1:length(imagefiles_ch1)
    
    currentfile = imagefiles_ch1(ii).name;
    filenames{ii} = erase(currentfile, '.tiff');
    currentimg_ch1 = imread(currentfile,1);
    % if using background correction:
    %images_ch1{ii} = wnCorr(currentimg_ch1);
    images_ch1{ii} = currentimg_ch1; 

end

for ii = 1:length(imagefiles_ch2)
    
    currentfile = imagefiles_ch2(ii).name;     
    currentimg_ch2 = imread(currentfile,1);
    images_ch2{ii} = currentimg_ch2;
    
end
    
ch1_parameters = {};
ch2_parameters = {};
iccs_parameters = {}; 

% analyzing pairs of channels in every image

for ii = 1:length(images_ch1)
    
    ch1 = images_ch1{ii};
    ch2 = images_ch2{ii};
    imgseries = cat(3, ch1, ch2);
    
    % calculating autocorrelation functions    
    ICS2DCorr_ch1 = corrfunc(imgseries(:,:,1));
    ICS2DCorr_ch2 = corrfunc(imgseries(:,:,2));
    
    % cropping and plotting each autocorrelation function 
    ICS2DCorrCrop_ch1 = autocrop(ICS2DCorr_ch1,12);
    figure(3);
    f = surf(ICS2DCorrCrop_ch1(:,:,1));
    axis tight
    colormap(jet)
    xlabel('\eta', 'FontSize', 12)
    ylabel('\xi', 'FontSize', 12)
    zlabel('r(\xi, \eta)', 'FontSize', 12)
    title('ICS2DCrop - channel 1')
    figname = strcat('ICS2DCrop_ch1',filenames{ii});
    saveas(f,figname,'png');
    
    ICS2DCorrCrop_ch2 = autocrop(ICS2DCorr_ch2,12);
    figure(4);
    g = surf(ICS2DCorrCrop_ch2(:,:,1));
    axis tight
    colormap(jet)
    xlabel('\eta', 'FontSize', 12)
    ylabel('\xi', 'FontSize', 12)
    zlabel('r(\xi, \eta)', 'FontSize', 12)
    title('ICS2DCrop - channel 2')
    figname = strcat('ICS2DCrop_ch2',filenames{ii});
    saveas(g,figname,'png');
    
    % storing each autocorrelation function
    ICS_ch1_series{ii} = ICS2DCorrCrop_ch1;
    ICS_ch2_series{ii} = ICS2DCorrCrop_ch2;
    
    % computing and cropping the cross-correlation function for each pair
    ICCS = corrfunc_cross(imgseries(:,:,1),imgseries(:,:,2));
    ICCS_crop = autocrop(ICCS,12);
    ICCS_series{ii} = ICCS_crop; 
    
    %plotting and saving the cross-correlation function
    h = surf(ICCS_crop(:,:,1));
    axis tight
    colormap(jet)
    xlabel('\eta', 'FontSize', 12)
    ylabel('\xi', 'FontSize', 12)
    zlabel('r(\xi, \eta)', 'FontSize', 12)
    set(h, 'LineStyle', 'none')
    title('ICCS2DCrop')
    figname = strcat('ICCS_crop',filenames{ii});
    saveas(h,figname,'png');
    
    % fitting the cross-correlation and two autocorrelation functions for
    % each image
    a = gaussfit(ICS2DCorrCrop_ch1, '2d', 0.066, 'y');
    ch1_parameters{ii} = a;
    plotgaussfit_bl(a(1,1:6), ICS2DCorrCrop_ch1(:,:,1), 0.066, 'y')
    figname = strcat('ICS2DCrop_ch1_fit',filenames{ii});
    savefig(figname);

    b = gaussfit(ICS2DCorrCrop_ch2, '2d', 0.066, 'y');
    ch2_parameters{ii} = b;
    plotgaussfit_bl(b(1,1:6), ICS2DCorrCrop_ch2(:,:,1), 0.066, 'y');
    figname = strcat('ICS2DCrop_ch2_fit',filenames{ii});
    savefig(figname);

    c = gaussfit(ICCS_crop, '2d', 0.066, 'y');
    iccs_parameters{ii} = c;
    plotgaussfit_bl(c(1,1:6), ICCS_crop(:,:,1), 0.066, 'y');
    figname = strcat('ICCS2DCrop_fit',filenames{ii});
    savefig(figname);

end

% calculating density of particles and the beam area from the results of
% Gaussian fitting

ch1_channel_fits = {};

for i = 1:length(ch1_parameters)
    params = cell2mat(ch1_parameters(i));
    particlesPerBeamArea_ch1 = 1/params(1);
    beamArea_ch1 = pi*params(2)*params(3);
    density_ch1 = particlesPerBeamArea_ch1/beamArea_ch1;
    ch1_channel_fits{i} = [particlesPerBeamArea_ch1 beamArea_ch1 density_ch1];
end

ch2_channel_fits = {};

for i = 1:length(ch2_parameters)
    params = cell2mat(ch2_parameters(i));
    particlesPerBeamArea_ch2 = 1/params(1);
    beamArea_ch2 = pi*params(2)*params(3);
    density_ch2 = particlesPerBeamArea_ch2/beamArea_ch2;
    ch2_channel_fits{i} = [particlesPerBeamArea_ch2 beamArea_ch2 density_ch2];
end

iccs_fits = {};

for i = 1:length(iccs_parameters)
    iccs_params = cell2mat(iccs_parameters(i));
    ch1_params = cell2mat(ch1_channel_fits(i));
    ch2_params = cell2mat(ch2_channel_fits(i));
    particlesPerBeamArea_iccs = (iccs_params(1)/((1/ch1_params(1))*(1/ch2_params(1))));
    beamArea_iccs = pi*iccs_params(2)*iccs_params(3);
    density_iccs = particlesPerBeamArea_iccs/beamArea_iccs;
    iccs_fits{i} = [particlesPerBeamArea_iccs beamArea_iccs density_iccs];    
end

% calculating colocalized fractions

M = [];

for i = 1:length(ch1_channel_fits)
    
    a = cell2mat(ch1_channel_fits(i));
    b = cell2mat(ch2_channel_fits(i));
    c = cell2mat(iccs_fits(i));
     
    M1 = c(:,1) / a(:,1);
    M2 = c(:,1) / b(:,1);
  
    M(:,:,i) = [a b c M1 M2]; 
    
end

M1 = M(:,10,:);
M2 = M(:,11,:);
    
%final plot

mean_M1 = mean(M1);
std_M1 = std(M1);
mean_M2 = mean(M2);
std_M2 = std(M2);

int = [mean_M1 mean_M2];
labels = {'M1' 'M2'};
error = [std_M1 std_M2];

figure
hold on
bar(int, 'FaceColor', [0.75 0.75 0.75]);
set(gca, 'XTickLabel', labels, 'XTick',1:numel(labels), 'fontsize', 14)
ylabel('Colocalized fractions')
err = errorbar(int,error, 'LineStyle', 'none');  
err.Color = [0 0 0];                            
%err.LineStyle = 'none'; 
hold off
