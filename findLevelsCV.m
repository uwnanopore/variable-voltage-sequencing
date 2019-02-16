function read = findLevelsCV(read, varargin)



sensitivity = 5;
minlevellength = 3;

for ca = 1:length(varargin)
    switch upper(varargin{ca})
        case 'SENSITIVITY'
            sensitivity = varargin{ca+1};
        case 'MINLEVELLENGTH'
            minlevellength = varargin{ca+1};
            
    end
end


eventdata = read.current_data;

badpts = filterFlickers(eventdata);

orig_ix = 1:length(eventdata);
eventdata = eventdata(~badpts);
orig_ix = orig_ix(~badpts);


[transitions, features, ~, stiffnesses] = findLevels_CVsub(eventdata, 'sensitivity', sensitivity, 'minlevellength', minlevellength);
starts = transitions(1:end-1)+1;
ends = transitions(2:end);
starts = orig_ix(starts);
ends = orig_ix(ends);


read.x_uf = features/read.OS;



for cL = 1:length(stiffnesses)
    stiffnesses{cL}(1,1) =  stiffnesses{cL}(1,1)*read.OS^2;
    stiffnesses{cL}(end,end) =  stiffnesses{cL}(end,end)*read.OS^2;
end

read.k_uf = stiffnesses;
read.npts_uf = ends-starts+1;
read.transitions = transitions(1:end-1);


end