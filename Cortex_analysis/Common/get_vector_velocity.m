function vectorout= get_vector_velocity(vectorin,params)


tracesmooth = medfilt1(vectorin,params.medfiltorder);
gfilter = fspecial('gaussian',[50 1], params.gaussorder);
tracesmoothed = convn(tracesmooth,gfilter,'same');
vectorout = cat(1,zeros((params.difforder)-2,1),tracesmoothed(params.difforder:end,1) - tracesmoothed(1:end-params.difforder+1),1);

end