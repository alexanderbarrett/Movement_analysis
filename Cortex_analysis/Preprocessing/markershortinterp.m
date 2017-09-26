function [interp_vals,fake_frames] = markershortinterp(markervals,recovery_length,max_spike)
markerind = 1;
gaps = cat(1,find(sum(markervals,2) == 0),find(isnan(sum(markervals,2))));
interp_vals = markervals;

% find all gaps and label with unique numbers
gap_base = zeros(size(markervals));
gap_base(gaps) = 1;
%[gaplabels,labelnum] = bwlabel(gap_base);
cc = bwconncomp(gap_base);
labelnum = cc.NumObjects;

%loop over gaps and get lengths
%could probably use arrayfun
gap_lengths = zeros(1,labelnum);
gap_border = zeros(1,labelnum);
fprintf('number of gaps %f \n',labelnum)
fake_frames = [];



for kk = 1:labelnum
    gap_lengths(kk) = numel(cc.PixelIdxList{kk});
    if (gap_lengths(kk)<recovery_length)
        gapstart = cc.PixelIdxList{kk}(1);
        gapend = cc.PixelIdxList{kk}(end);
        
        if (kk == 1 && labelnum>1)
            %make the gap border the difference between the neighboring gaps
            gap_border(kk) = min(gapstart-1,...
                cc.PixelIdxList{kk+1}(1)-...
                gapend -1);
            
            
        elseif (kk>1 && kk<labelnum )
            gap_border(kk) = min( gapstart-...
                cc.PixelIdxList{kk-1}(end),...
                cc.PixelIdxList{kk+1}(1)-...
                gapend -1);
            
        elseif (kk == labelnum && labelnum >1)
            gap_border(kk) = min( gapstart-...
                cc.PixelIdxList{kk-1}(end)-1,...
                length(markervals)- gapend -1);
            
        else
            gapborder(kk) = 3
        end
        
        %% do the interpolation
        xq = (gapstart:gapend);
        if gapstart> 1 && gapend < size(markervals,1)-1
            for jkj = 1:size(markervals,2)
                if (gap_border(kk) >3)
                    x = cat(2,gapstart-3:gapstart-1,gapend+1:(gapend+3));
                    v = markervals(x,jkj);
                    vq = interp1(x,v,xq,'spline');
                    
                else
                    x = cat(2,gapstart-1,gapend+1);
                    v = markervals(x,jkj);
                    vq = interp1(x,v,xq,'linear');
                    
                end
                interp_vals(xq,jkj) = vq;
                
            end
            fake_frames  = cat(2,fake_frames ,xq);
        end
    
end
end

%
% intergaps = find(gap_lengths<recovery_length);
% fake_frames = [];
% %if there is space on each side use a spline
% for mm = intergaps
%      gapstart = find(gaplabels == mm,1,'first');
%        gapend = find(gaplabels == mm,1,'last');
%
%        xq = (gapstart:gapend);
%       if gapstart> 1 && gapend < length(markervals)-1
%                  for jkj = 1:3
%
%    if (gap_border(mm) >3)
%        x = cat(2,gapstart-3:gapstart-1,gapend+1:(gapend+3));
%        v = markervals(x,jkj);
%        vq = interp1(x,v,xq,'spline');
%
%    else
%         x = cat(2,gapstart-1,gapend+1);
%        v = markervals(x,jkj);
%        vq = interp1(x,v,xq,'linear');
%
%    end
%           interp_vals(xq,jkj) = vq;
%
%       end
%      fake_frames  = cat(2,fake_frames ,xq);
%       end
% end
%
%
% % if not use linear interpolation
%
% %spikediff = diff(interp_vals);
% %difftimes = find(spikediff>max_spike);
%
%


end