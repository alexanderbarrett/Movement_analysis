function [fullposture_annotation_struct,fulloutput_annotation_struct,subset_of_frames_annotated,annotation_filenumber] = load_handannotation_files(annotation_number,mocapmasterdirectory)



gapfill_number = 20;

%% can load in multiple filenames
%annotation_mocapname =
%'Y:\Jesse\Data\Motionanalysis_captures\Vicon8\20170822\Preprocessed\nolj_Recording_day7_caffeine1_nolj.mat';.
switch annotation_number
    case 1
        %% FIRST FILE
        ratname = 'Vicon8';
        conditionname = 'Vicon8_caff';
        conditionnumber = 6;
        
        annotation_folder = strcat(mocapmasterdirectory,filesep,'Vicon8\20170822\Preprocessed\');
        
        annotation_filenames = {'nolj_Recording_day7_caffeine1_nolj_handannotation_JDM.mat'};
        annotation_filenumber = {1};
    case 2
        %% SECOND GROUP OF FILES
        ratname = 'Vicon8';
        conditionname = 'Vicon8_caff';
        conditionnumber = 6;
        
        annotation_folder = 'Y:\Jesse\Data\Motionanalysis_captures\Vicon8\20170822\Preprocessed\';
        annotation_filenames = {'nolj_Recording_day7_caffeine1_nolj_handannotation_JDM.mat',...
            'nolj_Recording_day7_caffeine1_nolj_handannotation_AB_1to48000.mat',...
            'nolj_Recording_day7_caffeine2_nolj_handannotation_AB_1to23800.mat',...
            'nolj_Recording_day7_caffeine2_nolj_handannotation_AB_23800to53323.mat',...
            'nolj_Recording_day7_caffeine2_nolj_handannotation_AB_53323to64603.mat'};
        
        annotation_filenumber = {1,1,2,2,2};
        
    case 3
        %% THIRD GROUP
        conditionname = 'Vicon8_prelesion';
        conditionnumber = 1;
        
        %annotation_folder = 'Y:\Jesse\Data\Motionanalysis_captures\Vicon8\20170820\Preprocessed\';
        %annotation_filenames = {'Recording_day5_overnight27_nolj_handannotation2_AB.mat'
        %annotation_filenumber = {23};
        
        annotation_folder = strcat(mocapmasterdirectory,filesep,'Vicon8\20170817\Preprocessed\');
        
        annotation_filenames = {'Recording_overnight_day28_nolj_handannotation_AB',...
            'Recording_overnight_day29_nolj_handannotation_AB',...
            'Recording_overnight_day211_nolj_handannotation_AB',...
            'Recording_overnight_day215_nolj_handannotation_AB',...
            'Recording_overnight_day216_nolj_handannotation_AB',...
            'Recording_overnight_day218_nolj_handannotation_AB',...
            'Recording_overnight_day220_nolj_handannotation_AB'...
            };
        annotation_filenumber = {1,2,4,8,9,11,13};
        
        
        %% fourth group
    case 4
        conditionname = 'JDM25_caff';
        conditionnumber = 5;
        
        annotation_folder = strcat(mocapmasterdirectory,filesep,'JDM25\20170921\Preprocessed\');
        
        annotation_filenames = {'nolj_Recording_day8_caff1_nolj_handannotation1_AB.mat',...
            'nolj_Recording_day8_caff1_nolj_handannotation2_AB',...
            'nolj_Recording_day8_caff2_nolj_handannotation_AB',...
            'nolj_Recording_day8_caff4_nolj_handannotation_AB'};
        annotation_filenumber = {1,1,2,4};
    case 5
        conditionname = 'JDM33_caff';
        conditionnumber = 1;
        
        annotation_folder = strcat(mocapmasterdirectory,filesep,'JDM33\20171130\Preprocessed\');
        
        annotation_filenames = {'nolj_Recording_day8_caff1_nolj_handannotation_AB.mat',...
            'nolj_Recording_day8_caff1_nolj_handannotation2_AB.mat'...
            };
        annotation_filenumber = {1,1};
    case 6
        %% NOT ACTUALLY ANNOTATED
        conditionname = 'JDM33_prelesion';
        conditionnumber = 1;
        
        annotation_folder = strcat(mocapmasterdirectory,filesep,'JDM33\20170921\Preprocessed\');
        
        annotation_filenames = {'nolj_Recording_day8_caff1_nolj_handannotation_AB.mat',...
            'nolj_Recording_day8_caff1_nolj_handannotation2_AB.mat'...
            };
        annotation_filenumber = {1,1};
    case 7
        conditionname = 'JDM25_prelesion';
        conditionnumber = 1;
        
        annotation_folder = strcat(mocapmasterdirectory,filesep,'JDM25\20170918\Preprocessed\');
        
        annotation_filenames = {'nolj_Recording_day5_overnight15_nolj_handannotation_AB.mat',...
            'nolj_Recording_day5_overnight5_nolj_handannotation_AB.mat'...
            };
        annotation_filenumber = {5,15};
    case 8
        conditionname = 'JDM25_postlesion';
        conditionnumber = 1;
        
        annotation_folder = strcat(mocapmasterdirectory,filesep,'JDM25\20171010\Preprocessed\');
        
        annotation_filenames = {'nolj_Recording_day27_overnight1_nolj_handannotation_AB.mat',...
            'nolj_Recording_day27_overnight3_nolj_handannotation_AB.mat'...
            };
        annotation_filenumber = {1,3};
    case 9
        conditionname = 'JDM33_postlesion';
        conditionnumber = 1;
        
        annotation_folder = strcat(mocapmasterdirectory,filesep,'JDM33\20180107\Preprocessed\');
        
        annotation_filenames = {'nolj_Recording_day45_overnight13_nolj_handannotation_AB.mat',...
            'nolj_Recording_day45_overnight14_nolj_handannotation_AB'...
            };
        annotation_filenumber = {13,14};
    case 10
        conditionname = 'JDM33_postlesion2';
        conditionnumber = 1;
        
        annotation_folder = strcat(mocapmasterdirectory,filesep,'JDM33\20180108\Preprocessed\');
        
        annotation_filenames = {'nolj_Recording_day46_overnight6_nolj_handannotation_AB.mat',...
            'nolj_Recording_day46_overnight7_nolj_handannotation_AB'...
            };
        annotation_filenumber = {6,7};
end


%% preprocess the annotation labels
annot_cell = cell(1,numel(annotation_filenames));
subset_of_frames_annotated = [];
agg_struct = struct([]);
agg_posture = struct([]);

for jj = 1:numel(annotation_filenames)%2:5
    output1 = load(strcat(annotation_folder,annotation_filenames{jj}));
    C = struct2cell(output1.output.GlobalBehavior);
    minval = min([C{:}]);
    maxval = max([C{:}]);
    
    for kk = fieldnames(output1.output.GlobalBehavior)'
        output1.output.GlobalBehavior.(kk{1}) = output1.output.GlobalBehavior.(kk{1})+(annotation_filenumber{jj}-1)*540000;
    end
    for kk = fieldnames(output1.output.Posture)'
        output1.output.Posture.(kk{1}) = output1.output.Posture.(kk{1})+(annotation_filenumber{jj}-1)*540000;
    end
    annot_cell{jj} = output1;
    
    
    subset_of_frames_annotated = cat(2,subset_of_frames_annotated ,(annotation_filenumber{jj}-1)*540000+(minval:maxval)); %depends on the file(s) loaded
    agg_struct = mergeStructs_JDM(agg_struct,output1.output.GlobalBehavior);
    agg_posture = mergeStructs_JDM(agg_posture,output1.output.Posture);
    
end
subset_of_frames_annotated = unique(subset_of_frames_annotated);

[pose_struct_out,globalbehavior_struct_out] = posture_behavior_gapfill(agg_posture,agg_struct);
%outputstruct_annotation = load(annotation_filename);


%% LOAD IN THE ANNOTATION
% fill the gaps in the structure
% in the future can write code to accomodate multiple annotation_files
[fulloutput_annotation_struct,indivbouts_annotation_struct] = fillannotationgaps_struct(globalbehavior_struct_out ,gapfill_number);
[fullposture_annotation_struct,indivbouts_posture_struct] = fillannotationgaps_struct(pose_struct_out ,gapfill_number);

%[fulloutput_annotation_struct,indivbouts_annotation_struct] = fillannotationgaps_struct(agg_struct ,gapfill_number);
%[fullposture_annotation_struct,indivbouts_posture_struct] = fillannotationgaps_struct(agg_posture ,gapfill_number);


