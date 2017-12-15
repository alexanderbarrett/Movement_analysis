function annotated_pose_struct = auto_annotate_pose(markerstruct_in,mocapstruct_in)

annotated_pose_struct = struct();
  %high rear
  annotated_pose_struct.high_rear = find(( markerstruct_in.SpineF(:,3)-markerstruct_in.SpineL(:,3))>50);
    annotated_pose_struct.very_high_rear = find(( markerstruct_in.SpineF(:,3)-markerstruct_in.SpineL(:,3))>80);

  %low rear -- shortens stance more and more
    annotated_pose_struct.low_rear = find(( markerstruct_in.HeadB(:,3)-markerstruct_in.SpineF(:,3))<-50);
  
  %l/r groom 
  annotated_pose_struct.RGroom = find(( markerstruct_in.SpineF(:,2)-markerstruct_in.SpineL(:,2))>20);
    annotated_pose_struct.LGroom = find(( markerstruct_in.SpineF(:,2)-markerstruct_in.SpineL(:,2))<-20);

    annotated_pose_struct.moving = mocapstruct_in.move_frames;
     % time_subset = mocapstruct_pre.move_frames;
 % time_subset2 = mocapstruct_post.move_frames;
  
[licktime,levertime] = gettasktime(mocapstruct_in);
annotated_pose_struct.leverind = find(levertime==1);

end