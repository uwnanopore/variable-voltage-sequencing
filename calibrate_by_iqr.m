function [m,b] = calibrate_by_iqr(data, caldata, doplots, databounds)

if nargin<4 && size(data,1)==3
    databounds = [0.05, 0.48;
                  0.10, 0.85;
                  0.20, 1.12];
elseif nargin<4 && size(data,1)==1
    databounds = [0.12, 0.85];
end

if nargin < 3
    doplots = 0;
end

 dataN = size(data, 2);
 calN = size(caldata, 2);
 calhist = histc(caldata', 0:0.02:1)'/calN;
 rawhist = diag(1./dataN)*histc(data', 0:0.02:1)';

bestm = ones(1, size(data, 1)); bestb = zeros(1, size(data, 1));


for ii = 1:size(data,1)
    bestm(ii) = interquartilerange(caldata(ii,:))/interquartilerange(data(ii,:));
    bestb(ii) = median(caldata(ii,:)) - bestm(ii)*median(data(ii,:));
end

besthist = diag(1./dataN)*histc((diag(bestm)*data + diag(bestb)*ones(size(data)))', 0:0.02:1)';

if doplots
    for ii = 1:size(data,1)
        figure, hold on
        plot(0:0.02:1, besthist(ii,:), 'linewidth', 2)
        plot(0:0.02:1, rawhist(ii,:), 'linewidth', 2)
        plot(0:0.02:1, calhist(ii,:), 'linewidth', 2)
        legend('calibrated', 'raw', 'calibrator')
        xlabel('Current (% IOS)')
        ylabel('normalized counts')
        title(['Feature ' num2str(ii) ', N = ' num2str(dataN) ', m = ' num2str(bestm(ii)) ', b = ' num2str(bestb(ii))])
    end
end

m = bestm; b = bestb;

end