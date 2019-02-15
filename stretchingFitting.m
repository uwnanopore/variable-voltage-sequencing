%%
% load the data
load('stretching_analysis_result.mat');

% define voltages
volts = [100:1:200];

% define reference voltage
V_ref = 180;

% calculate the cummulative shift
boot = zeros(size(analysis.boots, 2), size(analysis.boots{1}, 2));
for ii = 1:size(boot, 1)
    boot(ii, :) = analysis.boots{ii};
end
boot = cumsum([zeros(size(boot, 1), 1) boot], 2);
shift = mean(boot, 1);

% make 180 mV the reference voltage
shift = shift - shift(volts == V_ref);
for ii = 1:size(boot, 1)
    boot(ii, :) = boot(ii, :) - boot(ii, volts == V_ref);
end

% flip the sign to match the stretching math
shift = -1 .* shift;
boot = -1 .* boot;

% calculate the uncertainties based on the bootstrap
dshift = std(boot);

%%
% define N_ref: number of nucleotides in MspA @ 180 mV
N_ref = 12;

% define dN_ref
dN_ref = .5;

% calculate w_meas
w_meas = (N_ref + shift) / (N_ref);

% calculate dw_meas
dw_meas = sqrt( ((((-shift) ./ (N_ref .^ 2)) .* (dN_ref)) .^ 2) + (((1 ./ N_ref) .* (dshift)) .^ 2) );
weights = 1 ./ (dw_meas .^ 2);

% figure out betas
betas = volts ./ V_ref;

% define ex-fjc parameters
V_ref = V_ref ./ 1000; % convert mV to V
b = 1.5e-9; % persistence length in meters
S = 800e-12; % stretching modulus of ssDNA in Newtons
T = 22; % temperature in degC
kT = (1.38e-23) * (273 + T); % kT in Joules
lss = 0.69e-9; % nucleotide spacing of ssDNA
C = 6.24e18; % electron charges per coulomb
e = 1/C; % coulombs per electron charge

% define our model
w_model = @(beta, q) beta .* ( (((q .* V_ref .* b) - kT) .* (S + (q .* V_ref))) ./ (((beta .* q .* V_ref .* b) - kT) .* (S + (beta .* q .* V_ref))));

% choose which values of charge to test our fitting at
q_tests = (1.40:.0020:1.55) .* ((1e9) / (C)); % converts electron charge per nanometer into coulomb charge per meter

%%
% initialize storage for all the tested models
models = cell(1, length(q_tests));

% calculate all the models
for ii = 1:length(q_tests)
    models{ii} = w_model(betas, q_tests(ii));
end

% calculate the residuals on all the models
sse = models;
for ii = 1:length(sse)
    sse{ii} = sum((sse{ii} - w_meas) .^ 2);
end
sse = cell2mat(sse);

% calculate the error-weighted residuals on all models
wsse = models;
for ii = 1:length(wsse)
    wsse{ii} = sum(((wsse{ii}(dw_meas ~= 0) - w_meas(dw_meas ~= 0)) ./ (dw_meas(dw_meas ~= 0))) .^ 2);
end
wsse = cell2mat(wsse);

% find the index of the best model
[~, ix] = min(sse);
q_best = q_tests(ix);
model = models{ix};

% find the index of the best weighted model
[~, wix] = min(wsse);
q_wbest = q_tests(wix);
modelw = models{wix};

% calculate the best fit to the largest reasonable stretch
wup_meas = w_meas;
wup_meas(w_meas > 1) = w_meas(w_meas > 1) + dw_meas(w_meas > 1);
wup_meas(w_meas < 1) = w_meas(w_meas < 1) - dw_meas(w_meas < 1);

wsseup = models;
for ii = 1:length(wsseup)
    wsseup{ii} = sum(((wsseup{ii}(dw_meas ~= 0) - wup_meas(dw_meas ~= 0)) ./ (dw_meas(dw_meas ~= 0))) .^ 2);
end
wsseup = cell2mat(wsseup);
[~, wixup] = min(wsseup);
qup_wbest = q_tests(wixup);
modelwup = models{wixup};

% calculate the best fit to the smallest reasonable stretch
wdo_meas = w_meas;
wdo_meas(w_meas > 1) = w_meas(w_meas > 1) - dw_meas(w_meas > 1);
wdo_meas(w_meas < 1) = w_meas(w_meas < 1) + dw_meas(w_meas < 1);

wssedo = models;
for ii = 1:length(wssedo)
    wssedo{ii} = sum(((wssedo{ii}(dw_meas ~= 0) - wdo_meas(dw_meas ~= 0)) ./ (dw_meas(dw_meas ~= 0))) .^ 2);
end
wssedo = cell2mat(wssedo);
[~, wixdo] = min(wssedo);
qdo_wbest = q_tests(wixdo);
modlewdo = models{wixdo};

%%
stretching = struct;
stretching.source_data = analysis;
stretching.shift = shift;
stretching.dshift = dshift;
stretching.mV = [100:1:200];
stretching.w_meas = w_meas;
stretching.dw_meas = dw_meas;
stretching.betas = betas;
stretching.N_ref = N_ref;
stretching.dN_ref = dN_ref;
stretching.V_ref = V_ref;
stretching.b = b;
stretching.S = S;
stretching.T = T;
stretching.kT = kT;
stretching.C = C;
stretching.e = e;
stretching.lss = lss;
stretching.model = w_model;
stretching.charge_units = 'electrons per nanometer';
stretching.charge = q_wbest * C * 1e-9;
stretching.charge_high = qup_wbest * C * 1e-9;
stretching.charge_low = qdo_wbest * C * 1e-9;
stretching.best_fit = model;

save('stretching_fitting_result.mat', 'stretching');

