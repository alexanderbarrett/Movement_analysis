%load the right links and color scheme for the markerset
function [markercolor,links] =  day_loading_header(inputstring)




toy_dog = 0;
toydog = 0;
do_tsne = 0;
rat_hands = 0;
rat_elbows = 0;
rat_hands_offset = 0;
rat_nohands_noelbows = 0;
vicon3_singlehand = 0;

markercolor = {'b','b','r','r','r','m','g','g','k','k','k','k','k'};
links = {[1,2],[2,3],[3,4],[5,6],[4,5],[5,7],[5,8],[4,6],[3,6],[7,8]};

if (strcmp(inputstring,'vicon3_18marker'))
    %elbow is black
    markercolor = {'b','b','b',...
        'r','r','r',...
        'm','m',...
        'g','g',... %hips
        'y','y','y',... %L arm
        'w','w','w',... %R arm
        'g','g',...
        'b','b','k','k','k',};
    links = {[1,2],[2,3],[1,3],...
        [2,6],[1,6],[3,6],... %head to spine
        [6,4],[5,4],...
        [6,7],[7,8],[4,8],[4,7],[5,7],...
        [5,9],[5,10],...
        [11,12],[6,13],[6,14],[11,13],[12,13],...
        [14,15],[14,16],[15,16],...
        [9,18],[10,17]};
end



if (strcmp(inputstring,'vicon8_20marker'))
    %elbow is black
    markercolor = {'b','b','b',...
        'r','r','r',...
        'm','m',...
        'c','g',... %hips
        'y','y','y',... %L arm
        'w','w','w',... %R arm
        'g','c',...
        'c','g','k','k','k',};
    links = {[1,2],[2,3],[1,3],...
        [2,4],[1,4],[3,4],... %head to spine
        [4,5],[5,6],...
        [4,7],[7,8],[5,8],[5,7],[6,7],...
        [6,9],[6,10],...
        [11,12],[4,13],[4,14],[11,13],[12,13],...
        [14,15],[14,16],[15,16],...
        [9,18],[10,17],...%knees
        [18,19],[17,20]};
end



if (strcmp(inputstring,'vicon3_singlehand'))
    markercolor = {'b','b','b','r','r','r','m','m','g','g','y','y','k','k','k'};
    links = {[1,2],[2,3],[1,3],...
        [3,6],[6,5],[5,4],...
        [6,7],[7,8],[6,8],[4,8],[4,7],[5,7],[5,8],...
        [5,9],[5,10],...
        [6,11],[11,12]};
end

if (strcmp(inputstring,'rat_hands'))
    markercolor = {'b','b','r','r','r','m','g','g','y','y','y','y'};
    links = {[1,2],[2,3],[3,4],[4,6],[5,7],[5,6],[4,5],[5,7],[5,8],[4,6],[3,6],[7,8],[3,9],[3,11],[9,10],[11,12]};
end

if (strcmp(inputstring,'rat_hands_offset'))
    markercolor = {'b','b','r','m','r','m','r','g','g','y','y','y','y'};
    links = {[1,2],[2,3],[3,4],[4,6],[5,7],[8,9],[5,6],[3,5],[4,5],[6,7],[7,8],[7,9],[3,11],[3,10],[10,12],[11,13]};
end


if (strcmp(inputstring,'rat_elbows'))
    markercolor = {'b','b','r','m','r','m','r','g','g','y','y','k','k'};
    links = {[1,2],[2,3],[3,5],[5,7],[3,4],[4,5],[5,6],[6,7],[7,8],[7,9],[8,9],...
        [3,10],[3,11]};
end


if (rat_nohands_noelbows)
    markercolor = {'b','b','r','r','r','m','g','g','k','k','k','k','k'};
    links = {[1,2],[2,3],[3,4],[5,6],[4,5],[5,7],[5,8],[4,6],[3,6],[7,8]};
end

if (toy_dog)
    markercolor = {'b','b','r','m','r','g','g','y','y'};
    links = {[1,2],[2,3],[3,5],[3,4],[5,4],[3,6],[6,7],[3,9],[8,9]};
end
