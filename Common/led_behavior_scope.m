%check the frame differences
load('on_off_scope.txt')
load('on_off_behavior.txt')

max_scope = max(on_off_scope(:,2));
min_scope = min(on_off_scope(:,1));

max_behavior = max(on_off_behavior(:,2));
min_behavior = min(on_off_behavior(:,1));

ind_scope = [];
for i = 1:size(on_off_scope,1)
   ind_scope = cat(2,ind_scope,(on_off_scope(i,1):on_off_scope(i,2)));
end

ind_behavior = [];
for i = 1:size(on_off_behavior,1)
   ind_behavior = cat(2,ind_behavior,(on_off_behavior(i,1):on_off_behavior(i,2)));
end


edited_scope = setxor(min_scope:max_scope,ind_scope);
edited_behavior = setxor(min_behavior:max_behavior,ind_behavior);

dec_scope = unique(round(edited_scope./4));
dec_behavior = unique(round(edited_behavior./3));

min_length = min(length(dec_scope),length(dec_behavior));

difference = zeros(1,min_length);
for i=1:min_length
    difference(i) = dec_scope(i)-dec_behavior(i);
end

figure(2)
plot(difference)