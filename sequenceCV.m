function sequencedread = sequenceCV(read)

p_step = 0.93; p_skip = 0.07; p_ext = 0.2; p_bad = 0.01;
scorecutoff = -3;

p =[ p_step, p_skip*p_ext.^(0:10)];

map = load('pore_model_6mer_constant_voltage.mat');
map = map.map;
transition_info = load('transition_info_hel308_6mer.mat');
transition_info = transition_info.transition_info;


m = read.calibration_scale;
b = read.calibration_offset;
x=read.x_tf(1,read.startlevel:read.endlevel)*m + b;
K=read.k_tf(read.startlevel:read.endlevel);

for ii = 1:length(K)
    C = inv(K{ii});
    K{ii} = 1/C(1,1);
end

P=repmat(reshape(p,[1,1,12]), length(x), 4);
Pbad = p_bad*ones(1, length(x));
sequencedread = read;
sequencedread.sequence = calculateSequenceVV(x, K, P, Pbad, 'map', map, 'transitioninfo', transition_info, 'scorecutoff', scorecutoff);

       
end
