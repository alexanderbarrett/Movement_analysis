function checkrestmove(mocapstruct)
figure(33)
plot(mocapstruct.markers_preproc.HeadF)
test = zeros(1,size(mocapstruct.markers_preproc.HeadF,1));
test(mocapstruct.rest_frames)=1;
hold on
bar(200.*test)
hold off
end