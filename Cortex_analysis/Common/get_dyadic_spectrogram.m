function [dyad_out,fr,time_clustering] = get_dyadic_spectrogram(agg_features,opts)

%% set the number of voices/octaves for the GMM
numVoices = 10;
numOctave = 5;

%% set wavelet properties
wname = 'gaus1';
dt = 0.01;
f0 = centfrq(wname);

%% set the frequency for the scales -- these are used later on tog enerate the 
minfreq = 0.5;
minscale = f0./(minfreq*dt);
s0 = minscale*dt;

a0 = 2^(1/numVoices);
scales = s0*a0.^(0:numOctave*numVoices);

%% get the frequencies from the wavelet -- this is what you ultimately want
Freq = scal2frq(scales,wname,dt);
    freq_rev = fliplr(Freq);

%% get initial
mm=1;
        freq_subset = freq_rev(1+(mm-1)*numVoices:(mm)*numVoices);
        freq_delta_here = (max(freq_subset)-min(freq_subset))./numVoices;
        freq_range_here = min(freq_subset):freq_delta_here:max(freq_subset);
 [~,fr_temp,time_clustering,pc_spectrograms_temp] = spectrogram(agg_features(1,:),opts.clustering_window,...
            opts.clustering_overlap,freq_range_here,opts.fps);

agg_spectrograms = zeros(size(agg_features,1)*numOctave*size(pc_spectrograms_temp,1),size(pc_spectrograms_temp,2));

row_counter = 0;
row_delta = size(pc_spectrograms_temp,1);

for k=1:size(agg_features,1)
    
   %% get a multiresolution spectrogram
    fprintf('starting multiresolution spectrogram for feature %f \n',k);
    fr = [];
    for mm =1:numOctave
        row_counter = row_counter+1;
        freq_subset = freq_rev(1+(mm-1)*numVoices:(mm)*numVoices);
        freq_delta_here = (max(freq_subset)-min(freq_subset))./numVoices;
        freq_range_here = min(freq_subset):freq_delta_here:max(freq_subset);
        %      freq_range_here
        if sum(agg_features(k,:),2) ~=0
        [~,fr_temp,time_clustering,pc_spectrograms_temp] = spectrogram(agg_features(k,:),opts.clustering_window,...
            opts.clustering_overlap,freq_range_here,opts.fps);
        agg_spectrograms(((row_counter-1)*row_delta+1):(row_counter*row_delta),:) = pc_spectrograms_temp;
        else
                    agg_spectrograms(((row_counter-1)*row_delta+1):(row_counter*row_delta),:) = 0;

        end
       % pc_spectrograms{k} =   cat(1,pc_spectrograms{k},pc_spectrograms_temp);
        fr = cat(1,fr,fr_temp);
    end
   
end
num_fr = numel(fr);
%agg_spectrograms = cell2mat(pc_spectrograms'); %second dimension is time base

%% normalize the spectrograms
agg_spectrograms = log(agg_spectrograms);
agg_spectrograms(isinf(agg_spectrograms)) = -20;
agg_spectrograms(isnan(agg_spectrograms)) = -20;
dyad_out = agg_spectrograms;
end