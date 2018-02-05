%interspikeinterval check

day=2;

        raw_rasters_matrix = zeros(experiment_data{1}.length_trace{day},...
            length(experiment_data{1}.day_id{day}));%mean(experiment_data{day}.raw_rasters_matrix
        raw_raw_rasters_matrix = zeros(size(raw_rasters_matrix));
        raw_traces_matrix = zeros(size(raw_rasters_matrix));
        
        numICs = size(raw_rasters_matrix,2);
        length_trace =size(raw_rasters_matrix,1);
        isi_dist = cell(1,numICs);
        agg_isi = [];
        isi_under = zeros(1,numICs);
        
        for j=1:length(experiment_data{1}.day_id{day} )
            cellid = experiment_data{1}.day_id{day}(j);
            raw_rasters_matrix(:,j) = filter(ones(1,5),1,experiment_data{cellid}.events{day});
            raw_raw_rasters_matrix(:,j) = filter(ones(1,1),1,experiment_data{cellid}.events{day});
            raw_traces_matrix(:,j) = experiment_data{cellid}.trace{day};
            
            if (sum(raw_raw_rasters_matrix(:,j)) > 1)
            spiketimes = find(raw_raw_rasters_matrix(:,j) == 1);
            spikecirc = circshift(spiketimes,1);
            spikecirc(1) = 0;
            spikediff = spiketimes-spikecirc;
            spikediff = spikediff(2:end);
            spikediff(spikediff < 5) = [];
            isi_dist{j} = spikediff;
            agg_isi = [agg_isi spikediff'];
            isi_under(j) = numel(spikediff < 300);
            end
        end
        raw_traces_matrix(raw_traces_matrix < 0) = 0;
        
       [~,ind] = sort(isi_under,'descend');
        ind(1:50)
        
        figure(1)
        %hist(agg_isi,3000)
        
          result = 100;
    while (result > 0)      
prompt = 'Enter IC number to Examine [0 for exit] ';
result = input(prompt);
        if (result >0)
            figure(2)
            subplot(2,1,1)
            hist(isi_dist{result},10)
        subplot(2,1,2)
        plot(squeeze(raw_traces_matrix(:,result)))
        else 
            return
        end
    end