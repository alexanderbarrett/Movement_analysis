%htrfile = 'E:\Bence\Data\Motionanalysis_captures\Vicon3\20170721\skeletoncaptures\nolj_Vicon3_recording_amphetamine_vid2.htr';
htrfile = 'E:\Bence\Data\Motionanalysis_captures\Vicon3\20170721\Generated_HTR_files\nolj_Vicon3_recording_amphetamine_vid2_skel.htr';

skeleton = htrReadFile(htrfile);


skeleton.tree;
limbnames = fieldnames(skeleton);
limbstart = 15;
limbstop = numel(limbnames);
num_limbs = limbstop-limbstart;

rotation_matricies = cell(1,num_limbs);

limb_start = zeros(3,num_limbs);
limb_stop = zeros(3,num_limbs);
nframes = numel(skeleton.(limbnames{limbstart}).Tx);

limbposition = struct();


%% obtain rotation and translationmatricies
for kk = limbstart:limbstop;

 transmatrix =   reshape(cat(2,skeleton.(limbnames{kk}).Tx, skeleton.(limbnames{kk}).Ty, skeleton.(limbnames{kk}).Tz)',3,1,nframes);
 rotmatrix = eul2rotm(cat(2,deg2rad(skeleton.(limbnames{kk}).Rz),...
     deg2rad(skeleton.(limbnames{kk}).Ry),...
     deg2rad(skeleton.(limbnames{kk}).Rx )), skeleton.EulerRotationOrder);
 
 limbvec = zeros(3,1,nframes);
 limbvec(2,1,:) = skeleton.(limbnames{kk}).SF*1;
 
 
  limbposition.(limbnames{kk}) = cat(2,transmatrix,transmatrix+ mtimesx(rotmatrix,limbvec));
    
  
end

figure(23)


for mm = 1:10:1000
    for kk = limbstart:limbstop
   plot3( [squeeze(limbposition.(limbnames{kk})(1,1,mm)) ...
                        squeeze(limbposition.(limbnames{kk})(1,2,mm)) ],...
                    [squeeze(limbposition.(limbnames{kk})(2,1,mm)) ...
                        squeeze(limbposition.(limbnames{kk})(2,2,mm)) ],...
                       [squeeze(limbposition.(limbnames{kk})(3,1,mm)) ...
                        squeeze(limbposition.(limbnames{kk})(3,2,mm)) ],'linewidth',3)
                    hold on
    end
     
    hold off
    
       zlim([0 250])
            xlim([-350 350])
            ylim([-350 350])
                        set(gca,'XTickLabels',[],'YTickLabels',[],'ZTickLabels',[])
    pause(0.05)
end

