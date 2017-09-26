   agg_mse_full = cat(4,agg_mse{:});
   
   for mm =1:3
   for ll = 1:4
figure(120+ll+10*mm)
plot(squeeze(agg_mse_full([1,2],ll,mm,:))')

   end
   end