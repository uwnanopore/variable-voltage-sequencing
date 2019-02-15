function newdata = downsampleinmatlab(data,rate)
if rate==1
    newdata = data;
    return;
end
numdspts = floor(length(data)/rate);
data = data(1:numdspts*rate);

newdata = mean(reshape(data,rate,numdspts),1);

end