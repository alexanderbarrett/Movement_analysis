%% struct dynamic field names

A = [{'walk' };
    {'still'  };
    {'rear'}];

clear s
for i = 1:length(A)
    a = A(i);
    s.(A{i}) = i;
end

%% figure out how to iterate through structs
walk1.num_annotated = 10;
walk1.num_predicted = 15;
walk1.num_correct = 12;

walk2.num_annotated = 0;
walk2.num_predicted = 5;
walk2.num_correct = 2;

walk3.num_annotated = 7;
walk3.num_predicted = 5;
walk3.num_correct = 8;

still1.num_annotated = 3;
still1.num_predicted = 4; 
still1.num_correct = 5;

still2.num_annotated = 30;
still2.num_predicted = 40; 
still2.num_correct = 50;

still3.num_annotated = 13;
still3.num_predicted = 14; 
still3.num_correct = 15;

rear1.num_annotated = 2;
rear1.num_predicted = 3;
rear1.num_correct = 4;

rear2.num_annotated = 20;
rear2.num_predicted = 30;
rear2.num_correct = 40;

rear3.num_annotated = 12;
rear3.num_predicted = 13;
rear3.num_correct = 14;

s1.walk = walk1;
s2.walk = walk2;
s3.walk = walk3;

s1.still = still1;
s2.still = still2;
s3.still = still3;

s1.rear = rear1;
s2.rear = rear2;
s3.rear = rear3;

% fields = fieldnames(all_behaviors);
% for i = 1:numel(fields)
%   fprintf(fields{i})
%   fprintf('\n')
%   display(all_behaviors.(fields{i})) 
% end
% final_struct.walk.num_annotated = 0;
% final_struct.walk.num_predicted = 0;
% final_struct.walk.num_correct = 0;
% final_struct.still.num_annotated = 0;
% final_struct.still.num_predicted = 0;
% final_struct.still.num_correct = 0;
% final_struct.rear.num_annotated = 0;
% final_struct.rear.num_predicted = 0;
% final_struct.rear.num_correct = 0;

observed_behaviors = {'walk', 'still', 'rear'} ;


all_behaviors(1) = s1;
all_behaviors(2) = s2;
all_behaviors(3) = s3;

fields = {'num_annotated', 'num_predicted', 'num_correct'};

for i = 1:numel(observed_behaviors)   
    for j = 1:numel(fields)
        final_struct.(observed_behaviors{i}).(fields{j}) = 0;        
        for k = 1:length(all_behaviors)
            final_struct.(observed_behaviors{i}).(fields{j}) = final_struct.(observed_behaviors{i}).(fields{j}) + all_behaviors(k).(observed_behaviors{i}).(fields{j});   
        end
    end    
    final_struct.(observed_behaviors{i}).num_fp = final_struct.(observed_behaviors{i}).num_predicted - final_struct.(observed_behaviors{i}).num_correct ;
    final_struct.(observed_behaviors{i}).num_fn = final_struct.(observed_behaviors{i}).num_annotated - final_struct.(observed_behaviors{i}).num_correct ;
    final_struct.(observed_behaviors{i}).mc_rate = 1 - (final_struct.(observed_behaviors{i}).num_correct / final_struct.(observed_behaviors{i}).num_annotated) ;
    final_struct.(observed_behaviors{i}).fp_rate = 1 - (final_struct.(observed_behaviors{i}).num_fp - final_struct.(observed_behaviors{i}).num_predicted) ;
    final_struct.(observed_behaviors{i}).fn_rate = 1 - (final_struct.(observed_behaviors{i}).num_fn - final_struct.(observed_behaviors{i}).num_annotated) ;
end

beh_folds(1) = s1;
beh_folds(2) = s2;
beh_folds(3) = s3;

% iterate through structs, creating the desired final struct 
for i = 2:numel(observed_behaviors)   
    for j = 1:numel(fields)
        all_behaviors.(observed_behaviors{i}).(fields{j}) = 0;        
        for k = 1:length(beh_folds)
            all_behaviors.(observed_behaviors{i}).(fields{j}) = all_behaviors.(observed_behaviors{i}).(fields{j}) + beh_folds(k).(observed_behaviors{i}).(fields{j});
        end
    end      
    all_behaviors.(observed_behaviors{i}).num_fp = all_behaviors.(observed_behaviors{i}).num_predicted - all_behaviors.(observed_behaviors{i}).num_correct_pred ;
    all_behaviors.(observed_behaviors{i}).num_fn = all_behaviors.(observed_behaviors{i}).num_annotated - all_behaviors.(observed_behaviors{i}).num_correct_pred ;
    all_behaviors.(observed_behaviors{i}).mc_rate = 1 - (all_behaviors.(observed_behaviors{i}).num_correct_pred / all_behaviors.(observed_behaviors{i}).num_annotated) ;
    all_behaviors.(observed_behaviors{i}).fp_rate = all_behaviors.(observed_behaviors{i}).num_fp / all_behaviors.(observed_behaviors{i}).num_predicted ;
    all_behaviors.(observed_behaviors{i}).fn_rate = all_behaviors.(observed_behaviors{i}).num_fn / all_behaviors.(observed_behaviors{i}).num_annotated ;
end
    
%% make a table
sample = rand(3,3);
rowNames = {'a','b','c'};
colNames = {'x','y','z'};
sTable = array2table(sample,'RowNames',rowNames,'VariableNames',colNames)

unique(labels);



sTable = array2table(cfmat,'RowNames',,'VariableNames',colNames);

%% make confusion matrix


% test
hl = [3,3,8,8,3,4,6,6];
tbl = [3,4,8,8,4,4,3,6];

cfmat = zeros(8);

for n = 1:8
    for i = 1:(length(hl))
        if hl(i) == n
            num = tbl(i);
            cfmat(n,num) = cfmat(n,num) + 1
        end
    end
end

for i = 1:length(cfmat)
    cfmat2(i,:) = cfmat(i,:) / sum(cfmat(i,:));
end 
 
% implementation
cfmat = zeros(length(fieldnames_beh) - 1);

for n = 1:(length(fieldnames_beh) - 1)
    for i = 1:(length(hand_labels))
        if hand_labels(i) == n
            num = treebagger_labels(i);
            cfmat(n,num) = cfmat(n,num) + 1
        end
    end
end

%% remove transition from labels and frames to see if it improves treebagger

frames = tot_annotated_frames_goodtracking;
labels = tot_annotation_output_goodtracking;

bad_num = find(strcmp(tot_fieldnames_beh,'FspineSwap')==1) - 1;
indices = labels ~= bad_num;
labels_no_bad = labels(indices);
frames_no_bad = frames(indices);
labels = labels_no_bad;
frames = frames_no_bad;

tot_annotated_frames_goodtracking = frames;
tot_annotation_output_goodtracking = labels;

%% plot error vs num observations

error = zeros(1,length(observed_behaviors));
num_obs = zeros(1,length(observed_behaviors));
error_nt = zeros(1,length(observed_behaviors));
num_obs_nt = zeros(1,length(observed_behaviors));

for i = 2:length(observed_behaviors)
    error(i) = all_beh.(observed_behaviors{i}).mc_rate;
    num_obs(i) = all_beh.(observed_behaviors{i}).num_annotated;
    error_nt(i) = all_beh_nt.(observed_behaviors{i}).mc_rate;
    num_obs_nt(i) = all_beh_nt.(observed_behaviors{i}).num_annotated;
end
     
scatter(num_obs,error)

%% cost matrix

behs = fieldnames(all_behaviors);
obs_beh = [] ; 

% remove empty behaviors
for i = 1:numel(behs)  
  if all_behaviors.(behs{i}).num_annotated ~= 0   
     obs_beh = [obs_beh, behs(i)]; 
  end
end

cost_mat = zeros(length(obs_beh));
classnames = zeros(length(obs_beh), 1) ; 
for i = 1:length(obs_beh)    
    beh_number = find(strcmp(fieldnames_beh,obs_beh(i))==1) - 1 ;    
    classnames(i) = beh_number ;
    prior = all_behaviors.(obs_beh{i}).prior_prob;
    
    % cost 1 
    % cost_mat(i,:) = 1 / all_behaviors.(obs_beh{i}).prior_prob ;    
    
    % cost 2     
    
    if prior < 0.01 
        cost_mat(i,:) = 1.07 ;
    elseif 0.001 <= prior && prior < 0.01 
        cost_mat(i,:) = 1.06 ;
    elseif 0.01 <= prior && prior < 0.03 
        cost_mat(i,:) = 1.05 ;
    elseif 0.03 <= prior && prior < 0.05 
        cost_mat(i,:) =1.04 ;
    elseif 0.05 <= prior && prior < 0.1 
        cost_mat(i,:) = 1.03 ;
    elseif 0.1 <= prior && prior < 0.15 
        cost_mat(i,:) = 1.02 ;
    elseif 0.15 <= prior && prior < 0.2 
        cost_mat(i,:) = 1.01 ;
    else 
        cost_mat(i,:) = 1 ;
    end      

        
    % cost 3 
%     cost_mat(i,:) = sqrt(1 / all_behaviors.(obs_beh{i}).prior_prob) ;
    
    % cost 4 
%     n = 1 / all_behaviors.(obs_beh{i}).prior_prob ; 
%     cost_mat(i,:) = log(1000 / n) ;    
%     cost_mat(i,:) = 1 ; 
    cost_mat(i,i) = 0 ;    
end
classnames = cellstr(int2str(classnames)) ; 

[candidateframes,score,Mdl_test] = findsimilarframes_mlfeatures_multifeatures_2(ML_features,labels,frames,frames_to_annotate,classnames,cost_mat);

%% oversample rare categories

v = [1 2 3 4 1];
u = repelem(v,[1; 0; 3; 4]);

labels1 = boolean(labels == 56);
%%%%%%%

num_to_rep_vec = zeros(1,length(frames));
ps_num = find(strcmp(fieldnames_beh,'prone_still')==1) - 1 ;
ps_num_annotated = length(find(labels == ps_num)) ; 
for i = 1:(length(observed_behaviors))
    beh_number = find(strcmp(fieldnames_beh,observed_behaviors(i))==1) - 1 ;
    num_annotated = length(find(labels == beh_number)) ;
    if num_annotated ~= 0 
        indices = boolean(labels == beh_number) ;
        indices = round(indices * ps_num_annotated / num_annotated) ;
        num_to_rep_vec = num_to_rep_vec + indices ;
    end
end
oversampled_labels = repelem(labels,num_to_rep_vec) ; 
oversampled_frames = repelem(frames,num_to_rep_vec) ; 

labels = oversampled_labels ; 
frames = oversampled_frames ; 

  
%     if prior < 0.01 
%         indices = (indices * 50) ;
%     elseif 0.001 <= prior && prior < 0.01 
%         indices = (indices * 6) ;
%     elseif 0.01 <= prior && prior < 0.03 
%         indices = (indices * 6) ;
%     elseif 0.03 <= prior && prior < 0.05 
%         indices = (indices * 4) ;
%     elseif 0.05 <= prior && prior < 0.1 
%         indices = (indices * 2) ;
%     elseif 0.1 <= prior && prior < 0.15 
%         indices = (indices * 50) ;
%     elseif 0.15 <= prior && prior < 0.2 
%         indices = (indices * 50) ;
%     else 
%         indices = (indices * 50) ;
%     end  
    

num_to_rep_vec = zeros(1,length(frames));
ps_num = find(strcmp(fieldnames_beh,'Walk')==1) - 1 ;
ps_num_annotated = length(find(labels == ps_num)) ; 
for i = 1:(length(observed_behaviors))
    beh_number = find(strcmp(fieldnames_beh,observed_behaviors(i))==1) - 1 ;
    num_annotated(i) = length(find(labels == beh_number)) 
    
end

%% Undersample 

behavior_to_undersample = 'Walk' ; 

bad_num = find(strcmp(fieldnames_beh,behavior_to_undersample)==1) - 1 ;
i = 1 ;
% indicator to remove instance of behavior
j = 0 ; 
while i < length(labels) 
    if labels(i) == bad_num
        % skip this instance (don't remove)
        if j == 0 
            j = 1 ;
            aa = find(labels((i+1):end) ~= bad_num) ; 
            i = i + aa(1) ;             
        % remove this instance 
        else
            j = 0 ; 
            aa = find(labels((i+1):end) ~= bad_num) ; 
            i2 = i + aa(1) ; 
            labels(i:i2) = [] ;
            frames(i:i2) = [] ;
            i = i + 1 ;
        end
    else 
        i = i + 1 ; 
    end    
end
            
%% get num_annotated vector

num_annotated_vec = [] ; 
for i = 1:length(observed_behaviors)
    beh_num = find(strcmp(fieldnames_beh,observed_behaviors(i))==1) - 1 ;
    num_ann = sum(labels == beh_num) ; 
    if num_ann ~= 0 
        num_annotated_vec = [num_annotated_vec; num_ann] ; 
    end
end

%% fix cfmat for excel output (remove NaN)

% for cfmat_percent
cfmat = cfmat_percent ; 
indices = all(~isnan(cfmat)) ;
cfmat = cfmat(:,indices) ; 
cfmat = cfmat(indices,:) ; 
fb = tot_fieldnames_beh(2:end) ; 
vert = fb(indices) ; 
horiz = vert' ; 
cfmat_percent_fixed = cfmat ; 

% cfmat_percent take 2
cfmat = cfmat_percent ; 
indices = all(~isnan(cfmat)) ;
cfmat = cfmat(:,indices) ; 
indices2 = zeros(size(cfmat,1),1) ;
for i = 1:length(indices2)
    if sum(cfmat(i,:)) == 0
        indices2(i) = 1; 
    end
end
indices3 = logical(indices + ~(indices2')) ; 
cfmat = cfmat(indices3,:) ;
fb = good_fieldnames_beh(2:end) ;
vert = fb(indices3) ;
horiz = fb(indices) ;
cfmat_percent_fixed = cfmat ; 

% cfmat_percent take 3 
cfmat = cfmat_percent ; 
indices = all(~isnan(cfmat)) ;
cfmat = cfmat(:,indices) ; 
indices2 = zeros(size(cfmat,1),1) ;
for i = 1:length(indices2)
    if sum(cfmat(i,:)) == 0
        indices2(i) = 1; 
    end
end
indices3 = logical(~(indices2')) ; 
cfmat = cfmat(indices3,:) ;
fb = new_fieldnames_beh(2:end) ;
vert = fb(indices3) ;
horiz = fb(indices) ;
cfmat_percent_fixed = cfmat ; 
    
% for cfmat_rawnum (run cfmat_percent first)
cfmat = cfmat_rawnum ; 
cfmat = cfmat(:,indices) ; 
cfmat = cfmat(indices,:) ; 
cfmat_rawnum_fixed = cfmat ; 

% cfmat_rawnum take 2 / 3
cfmat = cfmat_rawnum ; 
cfmat = cfmat(:,indices) ; 
cfmat = cfmat(indices3,:) ; 
cfmat_rawnum_fixed = cfmat ; 

% produce heatmaps
    figure(1)
    heatmap_percent = heatmap(horiz,vert,cfmat_percent_fixed,'XLabel','Actual Class','YLabel','Predicted Class','CellLabelFormat','%.2g') ; 
    title('Test Data: JDM25 Caff; Training Data: Vicon8 Caff; Feature Set: annotation_ii');
% title('Cross Validation (30% Held Out); Data: Vicon8 Caff; Feature Set: annotation_ii');

figure(2)
heatmap_rawnum = heatmap(horiz,vert,cfmat_rawnum_fixed,'ColorScaling','scaledcolumns', ...
    'XLabel','Actual Class','YLabel','Predicted Class') ; 

%% produce heatmaps (with adjusted fieldnames_beh)

% for cfmat_percent
cfmat = cfmat_percent ; 
indices = all(~isnan(cfmat)) ;
cfmat = cfmat(:,indices) ; 
cfmat = cfmat(indices,:) ; 
vert = tot_fieldnames_beh(indices) ; 
horiz = vert' ; 
cfmat_percent_fixed = cfmat ; 

% cfmat_percent take 2
cfmat = cfmat_percent ; 
indices = all(~isnan(cfmat)) ;
cfmat = cfmat(:,indices) ; 
indices2 = zeros(size(cfmat,1),1) ;
for i = 1:length(indices2)
    if sum(cfmat(i,:)) == 0
        indices2(i) = 1; 
    end
end
indices3 = logical(indices + ~(indices2')) ; 
cfmat = cfmat(indices3,:) ;
fb = good_fieldnames_beh(2:end) ;
vert = fb(indices3) ;
horiz = fb(indices) ;
cfmat_percent_fixed = cfmat ; 

% cfmat_percent take 3 
cfmat = cfmat_percent ; 
indices = all(~isnan(cfmat)) ;
cfmat = cfmat(:,indices) ; 
indices2 = zeros(size(cfmat,1),1) ;
for i = 1:length(indices2)
    if sum(cfmat(i,:)) == 0
        indices2(i) = 1; 
    end
end
indices3 = logical(~(indices2')) ; 
cfmat = cfmat(indices3,:) ;
vert = good_fieldnames_beh(indices3) ;
horiz = good_fieldnames_beh(indices) ;
cfmat_percent_fixed = cfmat ; 
    
% for cfmat_rawnum (run cfmat_percent first)
cfmat = cfmat_rawnum ; 
cfmat = cfmat(:,indices) ; 
cfmat = cfmat(indices,:) ; 
cfmat_rawnum_fixed = cfmat ; 

% cfmat_rawnum take 2 / 3
cfmat = cfmat_rawnum ; 
cfmat = cfmat(:,indices) ; 
cfmat = cfmat(indices3,:) ; 
cfmat_rawnum_fixed = cfmat ; 

% produce heatmaps
figure(1)
heatmap_percent = heatmap(horiz,vert,cfmat_percent_fixed,'XLabel','Actual Class','YLabel','Predicted Class','CellLabelFormat','%.1g') ; 
title('Test Data: JDM25 Caff; Training Data: Vicon8 Caff; Feature Set: annotation_ii');
% title('Cross Validation (30% Held Out); Data: Vicon8 Caff; Feature Set: annotation_ii');

figure(2)
heatmap_rawnum = heatmap(horiz,vert,cfmat_rawnum_fixed,'ColorScaling','scaledcolumns', ...
    'XLabel','Actual Class','YLabel','Predicted Class') ; 
 
%% check renumbering

jumps = [] ;
for i = 1:length(tot_annotated_frames_goodtracking) - 1
    if tot_annotated_frames_goodtracking(i) > tot_annotated_frames_goodtracking(i + 1)
        jumps = [jumps ; {{i, tot_annotated_frames_goodtracking(i)}, {i + 1, tot_annotated_frames_goodtracking(i+1)}}];
    end
end

jumps = [] ; 
for i = 1:length(tot_ML_features.framelist_true) - 1
    if tot_ML_features.framelist_true(i) > tot_ML_features.framelist_true(i + 1)
        jumps = [jumps ; {tot_ML_features.framelist_true(i), tot_ML_features.framelist_true(i+1)}];
    end
end

behs = fieldnames(all_behaviors) ;
test_num = [] ;
test_label = [] ; 
for i = 1:numel(behs)
  if all_behaviors.(behs{i}).num_annotated ~= 0   
      test_label =  [test_label ; behs(i)] ;
      num = find(strcmp(tot_fieldnames_beh,behs(i))==1) - 1 ;
      test_num = [ test_num ; num ] ; 
      
  end
end

behs = fieldnames(all_behaviors) ;
train_num = [] ;
train_label = [] ; 
for i = 1:numel(behs)
  if all_behaviors.(behs{i}).num_annotated ~= 0   
      train_label =  [train_label ; behs(i)] ;
      num = find(strcmp(tot_fieldnames_beh,behs(i))==1) - 1 ;
      train_num = [ train_num ; num ] ; 
      
  end
end

% check to see if the number and fieldnames reassignment worked
behs = fieldnames(all_behaviors) ;
nums = {'Behavior', 'Training' , 'Testing' } ; 
for i = 1:numel(behs)
  if all_behaviors.(behs{i}).num_annotated ~= 0   
      num = find(strcmp(good_fieldnames_beh,behs(i))==1) - 1 ;
      num2 = find(strcmp(new_fieldnames_beh,behs(i))==1) - 1 ;
      nums = [ nums ; {(behs{i}),num , num2} ] ;
  end
end


