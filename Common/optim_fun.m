function F = optim_fun(x,xdata)
chunk_size = 80000;
num_chunks = floor(numel(xdata)./chunk_size);
F = zeros(1,numel(xdata));
num_chunks
for kk = 1:num_chunks
    index = (num_chunks*kk+1):(num_chunks*(kk-1));
F(index) = x(kk).*sin(2*pi*x(kk+2*num_chunks).*xdata(index)+x(kk+num_chunks)); %
end
