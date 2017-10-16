function [descriptor_struct_1,mocapfilearray,mocapfilestruct] = get_mocap_files(rat,save_tag,mocapmasterdirectory)


switch rat
    case 'Vicon8'
mocapfilestruct = loadmocapfilestruct('Vicon8',mocapmasterdirectory);
    case 'JDM25'
mocapfilestruct = loadmocapfilestruct('JDM25',mocapmasterdirectory);
end

%% get the desired files
descriptor_struct_1 = struct();

switch save_tag
    case 'Vicon8_prelesion'
   descriptor_struct_1.day = 5;
descriptor_struct_1.tag = 'overnight1';
descriptor_struct_1.cond = 'PreLesion';

good_inds = find(cellfun(@numel,strfind(mocapfilestruct.(descriptor_struct_1.cond).mocapfiles{descriptor_struct_1.day},descriptor_struct_1.tag)));
mocapfilearray = mocapfilestruct.(descriptor_struct_1.cond).mocapfiles{descriptor_struct_1.day}(good_inds);
    case 'Vicon8_caff'
descriptor_struct_1.tag = 'caff';
descriptor_struct_1.cond = 'PreLesion';
   descriptor_struct_1.day = 6;

good_inds = find(cellfun(@numel,strfind(mocapfilestruct.(descriptor_struct_1.cond).mocapfiles{descriptor_struct_1.day},descriptor_struct_1.tag)));
mocapfilearray = mocapfilestruct.(descriptor_struct_1.cond).mocapfiles{descriptor_struct_1.day}(good_inds);
 
    case 'Vicon8_amph'
        descriptor_struct_1.tag = 'amph';
descriptor_struct_1.cond = 'PreLesion';
   descriptor_struct_1.day = 7;
good_inds = find(cellfun(@numel,strfind(mocapfilestruct.(descriptor_struct_1.cond).mocapfiles{descriptor_struct_1.day},descriptor_struct_1.tag)));
mocapfilearray = mocapfilestruct.(descriptor_struct_1.cond).mocapfiles{descriptor_struct_1.day}(good_inds); 

 case 'Vicon8_dlslesion_early'
       descriptor_struct_1.tag = 'overnight_two1';
descriptor_struct_1.cond = 'UniLesion';
   descriptor_struct_1.day = 2;
good_inds = find(cellfun(@numel,strfind(mocapfilestruct.(descriptor_struct_1.cond).mocapfiles{descriptor_struct_1.day},descriptor_struct_1.tag)));
mocapfilearray = mocapfilestruct.(descriptor_struct_1.cond).mocapfiles{descriptor_struct_1.day}(good_inds);


 case 'Vicon8_dlslesion_late'
    descriptor_struct_1.tag = 'overnight';
descriptor_struct_1.cond = 'UniLesion';
   descriptor_struct_1.day = 7;
good_inds = find(cellfun(@numel,strfind(mocapfilestruct.(descriptor_struct_1.cond).mocapfiles{descriptor_struct_1.day},descriptor_struct_1.tag)));
mocapfilearray = mocapfilestruct.(descriptor_struct_1.cond).mocapfiles{descriptor_struct_1.day}(good_inds(1:10));



      case 'Vicon8_veh'
              descriptor_struct_1.tag = 'veh';
descriptor_struct_1.cond = 'PreLesion';
   descriptor_struct_1.day = 5;
good_inds = find(cellfun(@numel,strfind(mocapfilestruct.(descriptor_struct_1.cond).mocapfiles{descriptor_struct_1.day},descriptor_struct_1.tag)));
mocapfilearray = mocapfilestruct.(descriptor_struct_1.cond).mocapfiles{descriptor_struct_1.day}(good_inds);
% 
           case 'JDM25pre'
descriptor_struct_1.tag = 'overnight1';
descriptor_struct_1.cond = 'PreLesion';
   descriptor_struct_1.day = 6;
   good_inds = find(cellfun(@numel,strfind(mocapfilestruct.(descriptor_struct_1.cond).mocapfiles{descriptor_struct_1.day},descriptor_struct_1.tag)));
mocapfilearray = mocapfilestruct.(descriptor_struct_1.cond).mocapfiles{descriptor_struct_1.day}(good_inds);

 case 'JDM25post'
descriptor_struct_1.tag = 'overnight1';
descriptor_struct_1.cond = 'UniLesion';
   descriptor_struct_1.day = 6;
   good_inds = find(cellfun(@numel,strfind(mocapfilestruct.(descriptor_struct_1.cond).mocapfiles{descriptor_struct_1.day},descriptor_struct_1.tag)));
mocapfilearray = mocapfilestruct.(descriptor_struct_1.cond).mocapfiles{descriptor_struct_1.day}(good_inds);


end

