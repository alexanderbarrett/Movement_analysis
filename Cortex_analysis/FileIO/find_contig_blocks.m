function contigoutput = find_contig_blocks(timeinput,size_thresh)

values = zeros(1,numel(timeinput));
values(timeinput) = 1;
conn_comp = bwconncomp(values);


contigoutput = [];

for mm = 1:conn_comp.NumObjects-1
    
    if (conn_comp.PixelIdxList{mm+1}(1)-conn_comp.PixelIdxList{mm}(end))<50
        
        values(conn_comp.PixelIdxList{mm}(end):(conn_comp.PixelIdxList{mm+1}(1)))=1;
    end
end


conn_comp = bwconncomp(values);


    for mm = 1:conn_comp.NumObjects

   if numel(conn_comp.PixelIdxList{mm}) > size_thresh
       contigoutput = cat(1,contigoutput,conn_comp.PixelIdxList{mm});
       
   end
    
end

end

