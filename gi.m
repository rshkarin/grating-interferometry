% Inspired by Dr. Venera Weinhardt's implementation

clear;
tic;

path_darks = '/home/rshkarin/gr-data2/results-test/dark-fields/';
path_flats = '/home/rshkarin/gr-data2/results-test/flat-fields2/';
path_projs = '/home/rshkarin/gr-data2/results-test/projections/';

grating_steps = 4;
n_flats = 4;
n_darks = 5;
n_projs = 4
roi = [1300 355 750 1340]; 

disp('Processing Dark images');

for k = 1:n_darks
    fname = sprintf('%sradio-df-%d.tif', path_darks, k-1);
    disp(sprintf('Reading Dark-Filed file: %s\n', fname));
    darks(:,:,k) = imcrop(double(imread(fname)./4), roi);
end

mean_darks = median(darks,3);

disp('Processing Flats images');

for k = 1:n_flats
    fname = sprintf('%sradio-ff-%d.tif', path_flats, k-1);
    disp(sprintf('Reading Flat-Field file: %s\n', fname));
    flats(:,:,k) = imcrop(double(imread(fname)./4),roi);
    flats(:,:,k) = filter_im3(flats(:,:,k));
    flats(:,:,k) = flats(:,:,k) - mean_darks(:,:);
end

mean_flats = mean(mean(flats(:,:,:),1),2);
fft_mean_flats = fft(mean_flats);
[z,x,y] = size(fft_mean_flats);
[in,ind] = max(fft_mean_flats(1:floor(y/2)));
n_periods = ind + 1; 

fprintf('Number of periods = %d\n',n_periods-1);

fft_flats = fft(flats(:,:,:),[],3);
phase_flats(:,:) = angle(fft_flats(:,:,n_periods));
absorb_flats(:,:) = abs(fft_flats(:,:,1));
max_absorb_flats(:,:) = abs(fft_flats(:,:,n_periods));

disp('Processing Projections images');

for k = 1:floor(n_projs/grating_steps)
    fprintf('Reconstructing projection = %d\n', k-1);
    
    for p = 1:grating_steps
        index = (k-1) * grating_steps + p;
        fname = sprintf('%sradio-pr-%d.tif', path_projs, index-1);
        projs(:,:,index) = imcrop(double(imread(fname)./4),roi);
        projs(:,:,index) = filter_im3(projs(:,:,index));
        projs(:,:,index) = projs(:,:,index) - mean_darks(:,:);
    end
    
    fft_projs = fft(projs(:,:,:),[],3);
    phase_projs(:,:) = angle(fft_projs(:,:,n_periods));
    absorb_projs(:,:) = abs(fft_projs(:,:,1));
    max_absorb_projs(:,:) = abs(fft_projs(:,:,n_periods));
    
    absorb_contrast(:,:,k) = absorb_projs ./ absorb_flats;
    phase_contrast = phase_projs - phase_flats;
    phase_contrast = phase_contrast - median(median(phase_contrast));
    phase_contrast(:,:,k) = phase_contrast-2*pi*floor((phase_contrast+pi)/(2*pi));
    dark_contrast(:,:,k) = filter_im3((max_absorb_projs.*absorb_flats)./(max_absorb_flats.*absorb_projs));
    
    dark_contrast(isnan(dark_contrast)) = 0;
    dark_contrast(isinf(dark_contrast)) = 0;
    
    absorb_contrast(isnan(absorb_contrast)) = 0;
    absorb_contrast(isinf(absorb_contrast)) = 0;
    
    phase_contrast(isnan(phase_contrast)) = 0;
    phase_contrast(isinf(phase_contrast)) = 0;
    
    imtool(absorb_contrast,[]);
    imtool(phase_contrast,[]);
    imtool(dark_contrast,[]);
end 

toc;
