% initialize storage
vv = struct;
vv.read_number = zeros(1, 0);
vv.section_number = zeros(1, 0);
vv.cut_number = zeros(1, 0);
vv.alignment = cell(1, 0);
vv.accuracy = zeros(1, 0);
% loop over all reads
for cR = 1:100
    for cS = 1:100
        try
            read = load(['variableVoltageReads/read_' num2str(cR) '_section_' num2str(cS) '.mat']);
            read = read.read;
        catch
            continue
        end
        % loop over all cut sites
        for cC = 1:numel(read.cutsite_starts)
            % store read, section, cut numbers
            vv.read_number(end + 1) = cR;
            vv.section_number(end + 1) = cS;
            vv.cut_number(end + 1) = cC;
            % grab out called alignment to reference sequence
            vv.alignment{end + 1} = read.cut_alignment{cC};
            % grab out overall accuracy
            vv.accuracy(end + 1) = read.cut_accuracy(cC);
        end
    end
end
% load and store in random calling results
randoms = load('variable_voltage_random_results.mat');
randoms = randoms.variable_voltage_random_results;
vv.random_accuracy_samples = randoms.random_accuracy_samples;
vv.random_accuracy = randoms.random_accuracy;
% generate a single big alignment for overall result statistics
vv.big_alignment = cell2mat(vv.alignment);
% option to save
save_on = true;
if save_on
    variable_voltage_sequencing_results = vv;
    save('variable_voltage_sequencing_results.mat', 'variable_voltage_sequencing_results', '-v7.3');
end

