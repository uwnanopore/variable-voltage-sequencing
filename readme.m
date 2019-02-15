%% PIPELINE TAKING RAW VARIABLE-VOLTAGE CURRENT DATA READY TO SEQUENCE:
% ----------------------------------------------------------------------------
% ----------------------------------------------------------------------------
% COMPLETE PIPELINE FUNCTION:
% ----------------------------------------------------------------------------
% processVV.m
%   inputs:
%       current_data
%           time series current data for the event to be processed. can
%           either be "raw" (pre flicker filtering) or "cleaned" (post
%           flicker filtering). all data in this deposition is "cleaned"
%           and should be run with the is_raw input set to false (as is the
%           default behaviour)
%       voltage_data
%           time series voltage data for the event to be processed.
%       is_raw
%           boolean specifying whether the input data is "raw" (true) or
%           "cleaned" (false), as discussed above. default behavior sets
%           is_raw to false.
%   outputs:
%       read
%           a struct containing results from each phase of data processing,
%           culminating in the feature vectors (and associated stiffnesses)
%           for each level in the event after all filters have been applied
%           (removal, recombination, reorder).
%   related files:
%       processing_parameters_vv.mat
%           contains a struct with the data processing parameters used in
%           our work for each step of the variable voltage data processing
%           procedure.
%   subroutines:
%       processing_parameters_vv.m
%           script to generate and save the the
%           processing_parameters_vv.mat file
%       filterFlickers.m
%           removes flicker cycles from variable voltage data (see
%           Supplemental Note 3)
%       findLevelsPC.m
%           finds enzyme steps in variable voltage data using prinipal
%           component basis functions (see Supplemental Note 4)
%       capacitanceCompensation.m
%           removes the capacitive current signal in our variable voltage
%           data caused by the bilayer capacitance (see Supplemental Note
%           5)
%       makeUnfilteredFeatures.m
%           extracts feature vectors from the capacitance compensation
%           results (see Supplemental Note 6)
%       removalFilter.m
%           removes "bad" levels from the event's set off levels (see
%           Supplemental Note 8.1)
%       recombinationFilter.m
%           recombines levels that were measured more than twice, either
%           due to enzyme missteps or over calling during level finding
%           (see Supplemental Note 8.2)
%       reorderFilter.m
%           reorders out-of-order levels left over after previous filtering
%           based on level-to-level continuity information (see
%           Supplemental Note 8.3)
% 
%           
% FLICKER FILTERING:
% ----------------------------------------------------------------------------
% filterFlickers.m
%   inputs: 
%       data: time series (uncompensated) current data from variable
%             voltage experiment, downsampled over full voltage cycles
%             (i.e. down to 200 Hz)
%   varargin:
%       pflicker: set score thresholds for point removal
%       windows: set window sizes for point removal (in number of points)
%   outputs:
%       bad_cycles: indices of voltage cycles to be removed due to flickers
% 
% 
% CHANGE POINT DETECTION:
% ----------------------------------------------------------------------------
% findLevelsPC.m
%   inputs:
%       data: raw variable-voltage current measurements
%       vdata: raw variable-voltage voltage measurements
%       T: period of the cycling frequency, in data points (default = 250)
%       basis_fns: an N_pts x N_basis_functions array with each column
%                  representing one basis function to model the data. 
%                  for default v-v data, should use PC(:, 1:5).
%   varargin:
%       sensitivity: multiplier setting sensitivity to find change points.
%                    lower finds more change points. default is 4
%       minlevellength: minimum number of points allowed to seperate change
%                       points. default is the period length T.
%   outputs:
%       transitions: the locations (in data points) where change points
%                    were detected in the current data
%   related files:
%       principal_components_for_level_finding.mat
%           contains the PC vectors used as basis functions for change
%           point detection
%       cpic_fits.mat
%           contains fits correcting for multiple-testing fallacy in change
%           point detection
% 
% 
% CAPACITANCE COMPENSATION: 
% ----------------------------------------------------------------------------
% capacitanceCompensation.m
%   inputs:
%       v_data: voltage data for variable-voltage event
%       i_data: current data for variable-voltage event
%       transitions: locations (in data points) of the change points
%                    detected by findLevelsPC.m
%       ac_period: period of voltage cycle in data points (def = 250)
%       low_voltage: minimum voltage to consider in cap comp (def = 95)
%       high_voltage: maximum voltage to consider in cap comp (def = 205)
%       sample_frequency: frequency of sampled data per second (def = 50000)
%   outputs:
%       c_data: compensated current data (time series)
%       levels: contains info on I-V for each cap comped enzyme step
%       up_levels: cap comped I-V info for only up-swing in voltage
%       down_levels: cap comped I-V info for only down-swing in voltage
%       levels_nocc: un-comped I-V for each complete enzyme step
%       up_levels_nocc: un-comped I-V for only up-swing voltage
%       down_levels_nocc: un-comped I-V for only down-swing voltage
%       phase: overall phase of the voltage signal for the event
%       phase_index: the assigned index in voltage cycle for each sampled
%                    time point
%       Xraw_3: estimate of the 3 PCA features based on each half-cycle
%               measurement
%       Craw_3: estimate of the 3 PCA feature covariance based on
%               half-cycle measurements
%       Kraw_3: estimate of the 3 PCA feature stiffness based on half-cycle
%               measurements
%       well_conditioned: whether or not the 3 PCA feature covariances are
%                         well conditioned (i.e. minimum 6 independent
%                         measurements)
%       pc: the principal components used to reduce the data
%       evaluation_voltage: voltage values at which current is evaluated
%                           going into the PCA reduction
%   related files:
%       equal_x_spacing.mat
%           locations in voltage and DNA position for interpolation of the
%           IV curve for each level yielding equally spaced measurements
%           along the DNA position
%       principal_components.mat
%           PC vectors used as the 3 dimensional basis for reduced feature
%           vectors in downstream analysis
%   subroutines:
%       fitPrep2.m
%           prepares x, y data for fitting with Matlab
%       filterSpikes.m
%           removes anomolous electrical spikes found during PCA reduction
% 
% 
% PCA CONDUCTANCE FEATURE EXTRACTION:
% ----------------------------------------------------------------------------
% makeUnfilteredFeatures.m
%   inputs:
%       levels: the "levels" output from the previous cap comp step,
%               containing the complete I-V information for each enzyme
%               step (post cap comp).
%   outputs:
%       features_uf: n_features (101) x n_levels array with the full
%                    (non-pca'd) feature vectors
%       stiffnesses_uf; a 1 x n_levels cell array, each containing the full
%                       inverse covariance matrix (diagonal) for the
%                       feature vector for each state
%       npts_uf: total number of data points in each level
%       evaluation_voltage: set of voltage values (in mV) where the
%                           features are extracted
%       evaluation_x: set of corresponding dna positions  (in nt)
%       feature_type: a string describing the nature of the extracted
%                     feature vectors
%   related files:
%       equal_x_spacing.mat
%           locations in voltage and DNA position for interpolation of the
%           IV curve for each level yielding equally spaced measurements
%           along the DNA position
%   
% 
% REMOVAL FILTER:
% ----------------------------------------------------------------------------
% removalFilter.m
%   inputs:
%       x_101: 101 x N_levels array of un-normalized, un-calibrated 101
%              dimensional feature vectors for event to be filtered
%   varargin:
%       threshold: sets probability threshold for level to be called as
%                  "bad". default is set to the logit(0), corresponding to
%                  the SVM boundary. user input can range from 0 (must have
%                  P(good) = 0 to be called "bad") to 1 (must have P(good)
%                  <= 1 to be called "bad"). Claser to 1, more levels will
%                  be called "bad". default corresponds to ~0.45.
%   outputs:
%       keep_ix: 1 x N_levels array of 1 (good levels to keep) and 0 (bad
%                levels to throw)
%       uf_to_f: 1 x N_levels array of alignment of uf levels to new f
%                levels
%       threshold: retures P(good) required to be below to call level "bad"
%   related files:
%       pore_model_6mer_variable_voltage.mat
%           6-mer pore model for variable voltage with the half-stepping
%           hel308
%       principal_components.mat
%           PC vectors used as the 3 dimensional basis for reduced feature
%           vectors in downstream analysis
%       svm_for_bad_level_filtering.mat
%           SVM for classifying states as "bad" levels for removal
%       logit_params_for_bad_level_filtering.mat
%           parameters for the logistic function convertind distance from
%           the bad level SVM decision boundary to probability of good/bad
%       logit_for_badlevel_filter.mat
%           anonymous function of the logit curve
%       removal_svm_training_examples.mat
%           training examples used to develop the removal filter SVM and
%           logit function
%   subroutines:
%       levelPredPipeline.m
%           lookup function using a pore model to predict the signal for a
%           given input DNA sequence
%       calibrate_by_iqr.m
%           calibrates two sets of data using a scale and offset to match
%           up their inter quartile ranges
%       buildRemovalSVM.m
%           script to construct the SVM and logit used in the removal
%           filter
%       trainRemovalClassifierQuadSVM.m
%           function conducting the actual SVM training for the removal SVM
% 
% 
% RECOMBINATION FILTER:
% ----------------------------------------------------------------------------
% recombinationFilter.m
%   inputs:
%       x_101: the un-calibrated, un-normalized, removal-filtered features
%              (non-pca) for the event
%       k_101: stiffnesses for the full features
%       x_3: pca-reduced feature vectors
%       k_3: stiffnesses for pca-reduced feature vectors
%       npts: number of raw data points contributing to each level
%   varargin:
%       enzymestepprobabilities: user input of priors on enzyme transition
%                                probabilities, will be folded into the
%                                contiguity-based SVM step predictions
%       principlecomponents: user input of de-noising princomps
%       calibrationmap: user input of kmer map to calibrate to for translation
%                       into right scaling for svm
%       referencesequence: sequence to calibrate measurements to
%       lookback: how many levels to look back in allowing alignment
%   outputs:
%       f_to_tf: alignment matching removal filtered levels to
%                recombination filtered levels
%       x_101_tf: full feature vectors after reco-filter
%       k_101_tf: stiffnesses of full feature vectors
%       x_3_tf: pca-reduced feature vectors after reco-filter
%       k_3_tf: stiffnesses of pca-reduced feature vectors
%       c_3_tf: covariances of pca-reduced feature vecotrs
%       npts_tf: number of data points in each level
%       did_backstep_f: boolean whether each input level was observed to
%                       backstep
%       did_backstep_tf: boolean whether each output level was observed to
%                        backstep
%   related files:
%       principal_components.mat
%           PC vectors used as the 3 dimensional basis for reduced feature
%           vectors in downstream analysis
%       svm_for_recombination_filter.mat
%           SVM for classifying transition types for recombination filter
%       logit_params_for_recombination_filtering.mat
%           logit function parameters for assigning probabilities to
%           different step types based on distance from SVM boundaries
%       logit_for_backstep_filter.mat
%           anonymous function for logit
%       validation_for_recombination_filter.mat
%           results on hold-out validation set to find overall accuracy of
%           the recombination svm, and set regularization ceiling on
%           assigned probabilities
%       pore_model_6mer_variable_voltage.mat
%           6-mer pore model for variable voltage with the half-stepping
%           hel308
%       training_data_for_recombination_filter.mat
%           data set used to train the recombination SVM
%       check_levels_for_recombination_filter.mat
%           boolean of which levels have been vetted for inclusion in the
%           training set
%       use_levels_for_recombination_filter.mat
%           boolean of which levels have been vetted and chosen for
%           inclusion in the training set
%       phix.fasta
%           fasta file giving the phix-174 dna sequence
%   subroutines:
%       levelPredPipeline.m
%           lookup function using a pore model to predict the signal for a
%           given input DNA sequence
%       normalizerAC.m
%           function to remove the non-dna-dependent component from the
%           conductance signal, allowing recovery of the continuous curve
%           as a function of DNA position. the pore model exists in
%           "normalized" space, so data is normalized to calibrate it
%           against the pore model
%       calibrate_by_iqr.m
%           calibrates two sets of data using a scale and offset to match
%           up their inter quartile ranges
%       calculateStepProbabilitiesForRecombinationFilter.m
%           function using reco SVMs and logits to calculate transition 
%           probabilities based on input feature vectors
%       selfAlignAC.m
%           function actually carrying out self alignment
%       smatrixACBF.m
%           function calculating the similarity scores between all variable
%           voltage levels for self alignment. allows an additional penalty
%           to be taken on comparison while throughing out a "bad" feature
%       multivarSampleStats.m
%           combines two levels into one, averaging according to their 
%           respective numbers of data points. returns a new mean and
%           stiffness matrix
%       buildRecombinationSVM.m
%           script training set of SVMs and logits for recombination filter
%       fitprep.m
%           prepares x, y data for fitting with Matlab
%       trainClassifierQuadSVM.m
%           trains a SVM classifier with quadratic kernel to seperate input
%           labeled data
%       
% 
% REORDER FILTER:
% ----------------------------------------------------------------------------
% reorderFilter.m
%   inputs:
%       x_3_tf: pca-reduced feature vectors of the reco-filtered data
%       k_3_tf: stiffnesses of pca-reduced feature vectors of the
%               reco-filtered data
%       x_101_tf: full feature vectors of the reco-filtered data
%       k_101_tf: stiffnesses of full feature vectors
%       npts_tf: number of points in each reco-filtered level
%       did_backstep_tf: boolean of whether each reco-filtered level was
%                        observed to take a backstep
%   varargin:
%       maxiterations: sets how many itererations of reordering to allow
%       map: sets the kmer model to calibrate to
%       prior: sets a prior on relative probabilities of step/backstep/skip
%   outputs:
%       x_3_r: pca-reduced feature vectors of the reord-filtered data
%       k_3_r: stiffnesses for pca-reduced reord-filt feature vectors
%       x_101_r: full feature vectors of reord-filt data
%       k_101_r: stiffnesses of full feature vectors
%       npts_r: number of points in each reord-filt level
%       did_backstep_r: boolean whether each reord-filt level was observed
%                       to take a backstep
%   related files:
%       principal_components.mat
%           PC vectors used as the 3 dimensional basis for reduced feature
%           vectors in downstream analysis     
%       pore_model_6mer_variable_voltage.mat
%           6-mer pore model for variable voltage with the half-stepping
%           hel308
%       svm_for_recombination_filter.mat
%           SVM for classifying transition types for recombination filter
%           ** use same classifier for reorder as for recombination **
%       logit_for_backstep_filter.mat
%           anonymous function for logit
%           ** use same classifier for reorder as for recombination **
%       logit_params_for_recombination_filtering.mat
%           logit function parameters for assigning probabilities to
%           different step types based on distance from SVM boundaries
%           ** use same classifier for reorder as for recombination **
%       validation_for_recombination_filter.mat
%           results on hold-out validation set to find overall accuracy of
%           the recombination svm, and set regularization ceiling on
%           assigned probabilities
%           ** use same classifier for reorder as for recombination **
%   subroutines:
%       levelPredPipeline.m
%           lookup function using a pore model to predict the signal for a
%           given input DNA sequence
%       normalizerAC.m
%           function to remove the non-dna-dependent component from the
%           conductance signal, allowing recovery of the continuous curve
%           as a function of DNA position.
%       calibrate_by_iqr.m
%           calibrates two sets of data using a scale and offset to match
%           up their inter quartile ranges
%       calculateStepProbabilitiesForReorderFilter.m
%           function using reordering SVMs and logits to calculate transition 
%           probabilities for reordering filter based on input feature vectors
%       calculateReordering.m
%           function implementing the dynamic programming algorithm to find
%           most likely set of transitions linking the observed states
% 
% 
%% CODE FOR CONDUCTING VARIABLE VOLTAGE SEQUENCING
% ----------------------------------------------------------------------------
% ----------------------------------------------------------------------------
% WRAPPER TO CONDUCT VARIABLE VOLTAGE SEQUENCING
% ----------------------------------------------------------------------------
% sequenceVV.m
%   purpose:
%       Serves as a wrapper, taking a processed variable voltage read (all
%       filtering steps conducted, and calibration applied) and determines
%       the maximum likelihood sequence to have generated it, subject to
%       allowed state-to-state transitions and a pore model mapping
%       sequence to signal.
%   inputs:
%       read:
%           a fully processed variable voltage read (as found in
%           /variableVoltageReads/read_*_section_*.mat). must contain the
%           fields x101_r, x3_r, x3_cal, k3_cal, calibration_scale, and
%           calibration_offset
%   outputs:
%       sequence:
%           maximum likelihood sequence generating the observed variable
%           voltage levels, given in 5' -> 3' order
%   related files:
%       pore_model_6mer_variable_voltage.mat
%           6-mer pore model for variable voltage with the half-stepping
%           hel308
%       transition_info_hel308_6mer.mat
%           information detailing which kmer model states can transition to
%           which other kmer model states, and with which associated step
%           penalties, assuming a 6-mer pore model and hel308 half stepping
%   subroutines:
%       smartStepCounts.m
%           At each enzyme step, calculates the probabilities of different
%           transition types (step, skip1, skip2, ..., bad) based on the
%           variable-voltage continuity information, similar as to how has been
%           done in recombination and reordering filters. Outputs the
%           step/skip/bad probabilities to be fed into the sequencer.
%       calculateSequenceVV.m
%           uses an HMM to determine the maximum likelihood sequence generating
%           an observed set of variable voltage levels.
%   
% CALCULATE TRANSITION PENALTIES
% ----------------------------------------------------------------------------
% smartStepCounts.m
%   purpose:
%       At each enzyme step, calculates the probabilities of different
%       transition types (step, skip1, skip2, ..., bad) based on the
%       variable-voltage continuity information, similar as to how has been
%       done in recombination and reordering filters. Outputs the
%       step/skip/bad probabilities to be fed into the sequencer.
%   inputs:
%       x_101
%           (101 x n_levels) features array. should be unnormalized and
%           uncalibrated. typically final sequencing will use the _r
%           features (reordered product)
%       x_3
%           (3 x n_levels) features array of 3-pc reduced feature vectors
%           corresponding to the full feature vectors in x_101
%       calibrated: 
%           (boolean) tells whether passed features have yet been
%           calibrated
%       calibration
%           (1 x 2) array of [scale, offset] to use if passed
%           "calibrated" is true
%   outputs:
%       p
%           (n_levels x max_sequential_bad x 12) array giving transition
%           probabilities into each level. 2nd dimension gives gap sizes,
%           3rd dimension gives step sizes. max sequential bad is default at
%           3, meaning at most the sequencer is allowed to discard 3
%           consecutive measured levels as "bad". more allowed sequential
%           bads could potentially improve performance but at the cost of
%           run time.
%       p_bad
%           (n_levels x 1) array of bad probabilities of each level
%   varargin:
%       dumb_mode
%           turns on dumb_mode where we just use default step
%           probabilities at each
%       defaultpbad
%           allows user to provide the default p_bad value
%       defaultpskip
%           allows user to provide the default p_skip value
%       defaultpextension
%           allows user to provide the default p_extend
%       maxsequentialbad
%           provide the max consecutive bad levels allowed to be called
%       regularization
%           turn on (true) or off (false) the logit regularization by
%           verification calling rates
%       pfalloff
%           fall off probability rate for extending steps after we no
%           longer trust the svms (0 = max fall off, 1 = no fall off)
%       priortostep
%           prior probability on how often steps should occur relative to
%           skips. default is 0.5, in which case we trust exactly what
%           comes out of the svm ladder. a value above 0.5 biases towards
%           more steps
%       map
%           allows user to pass in a map to save on loading time
%   related files:
%       classifiers_for_smart_step_counts.mat
%           set of SVM + logit classifiers for determining smart step
%           counts from input level data
%       smart_step_counts_training_examples.mat
%           set of examples used to train the smart step counts SVM + logit
%           classifiers
%       svm_for_bad_level_filtering.mat
%           SVM for classifying states as "bad" levels for removal
%       logit_params_for_bad_level_filtering.mat
%           parameters for the logistic function convertind distance from
%           the bad level SVM decision boundary to probability of good/bad
%       validation_for_bad_level_filtering.mat
%           results on hold-out validation set to find overall accuracy of
%           the removal svm, and set regularization ceiling on
%           assigned probabilities
%       logit_for_backstep_filter.mat
%           anonymous function of the logit curve
%       pore_model_6mer_variable_voltage.mat
%           6-mer pore model for variable voltage with the half-stepping
%           hel308
%   subroutines:
%       buildSmartStepCountsSVM.m
%           script training set of SVMs and logits for smart step counts
%       trainClassifierQuadSVM.m
%           trains a SVM classifier with quadratic kernel to seperate input
%           labeled data
%       levelPredPipeline.m
%           lookup function using a pore model to predict the signal for a
%           given input DNA sequence
%       normalizerAC.m
%           function to remove the non-dna-dependent component from the
%           conductance signal, allowing recovery of the continuous curve
%           as a function of DNA position.
%       calibrate_by_iqr.m
%           calibrates two sets of data using a scale and offset to match
%           up their inter quartile ranges
% 
%
% CONDUCT HMM SEQUENCING
% ----------------------------------------------------------------------------
% calculateSequenceVV.m
%   purpose:
%       uses an HMM to determine the maximum likelihood sequence generating
%       an observed set of variable voltage levels.
%   inputs:
%       x
%           the d x n array of features (d = dimension of feature vectors,
%           n = number of levels in event)
%       K
%           the 1 x n cell array of the d x d stiffness matrices on each
%           feature vector
%       p
%           an n x (b + 1) x 3 array whose 3rd dim elements are [p(step),
%           p(skip), p(skip extend)]. for example, p(12, 3, 1) is the step
%           probability from the 12 - 3 = 9th level to the 12th level. the
%           2nd dimension size is set by the maximum number of sequential
%           bad levels allowed to be found by the sequencer; it is 1 less
%           than the dimension (equal to b in the above expression)
%       p_bad
%           an n x 1 array of bad level probabilities for each level
%   outputs:
%       sequence
%           the maximum likelihood sequence generating the observed levels,
%           in 5' -> 3' order
%       total_score
%           the total accumulated alignment score of the sequencing run
%           against the pore model
%       kmer_calls
%           called map state matched by each level
%       kmer_confidence
%           confidence in choice of each map kmer matched to
%       kmer_half_steps
%           the called step size between each pair of levels
%       level_was_good
%           boolean whether or not the level was "good" and included in the
%           sequencing or deemed "bad" and effectively ignored
%       source_levels
%           measured level most responsible for each individual base call
%   varargin:
%       backsteps
%           an n x 1 logical array determining whether or not each level
%           had been observed to have backstepped (during previous
%           filtering analysis)
%       map
%           user input of the pore model to be used in decoding
%       transitioninfo
%           user input of the structure of transition information (which
%           says which map states can transition to which other map states,
%           and with which associated penalties)
%       pindback
%           user input of the probability that a level is an
%           ATP-independent hel308 step given that we observed it to
%           backstep (default = 0.975)
%       pback
%           overall probability of a backstep (default = 0.025)
%       scorecutoff
%           how much worse a match score can be than the best scoring match
%           in a row for it to still be considered. this is a negative
%           number. making its magnitude larger results in a more
%           comprehensive search, but takes longer. for complete alignment
%           against the pore model, use -inf. default value is -10
%   related files:
%       transition_info_hel308_6mer.mat
%           information detailing which kmer model states can transition to
%           which other kmer model states, and with which associated step
%           penalties, assuming a 6-mer pore model and hel308 half stepping
%       pore_model_6mer_variable_voltage.mat
%           6-mer pore model for variable voltage with the half-stepping
%           hel308
%   subroutines:
%       smatrixACpar.m
%           function calculating the similarity scores between all variable
%           voltage levels and the pore model levels. parallelized using Matlab
%           parallel computing toolbox to improve speed
%       logdet.m
%           computes the log determinant of a symmetric, positive definite matrix A
%
%
%% CODE FOR GENERATING MAIN TEXT FIGURES
% ----------------------------------------------------------------------------
% ----------------------------------------------------------------------------
% FIGURE 1
% ----------------------------------------------------------------------------
% Figure1.m
%   related files:
%       figure1_jumps_data.mat
%           contains data for the level transitions in figure 1b
%       figure1_reduced_data.mat
%           contains the current data plotted in figure 1b
%       pore_model_6mer_variable_voltage.mat
%           6-mer pore model for variable voltage with the half-stepping
%           hel308
%       pore_model_6mer_constant_voltage.mat
%           6-mer pore model for constant voltage with the half-stepping
%           hel308
%       stretching_fitting_result.mat
%           results of DNA elongation analysis that are plotted in figure
%           1e
%       equal_x_spacing.mat
%           locations in voltage and DNA position for interpolation of the
%           IV curve for each level yielding equally spaced measurements
%           along the DNA position
%   subroutines:
%       downsampleinmatlab.m
%           vectorized function to downsample time series data
%       levelPredPipeline.m
%           lookup function using a pore model to predict the signal for a
%           given input DNA sequence (variable voltage)
%       levelPredPipelineDC.m
%           lookup function using a pore model to predict the signal for a
%           given input DNA sequence (constant voltage)
%       shadedErrorBar.m
%           plotting script to generate images as in panel 1e
%     
% 
% FIGURE 2
% ----------------------------------------------------------------------------
% Figure2.m
%   related files:
%       principal_components.mat
%           PC vectors used as the 3 dimensional basis for reduced feature
%           vectors in downstream analysis   
%       pore_model_6mer_variable_voltage.mat
%           6-mer pore model for variable voltage with the half-stepping
%           hel308
%       figure2_eventdata_1.mat
%           contains data for plotting panels (a) and (c)
%       figure2_eventdata_2.mat
%           contains data for plotting panels (b) and (d)
%       svm_for_bad_level_filtering.mat
%           SVM for classifying levels as good or bad
%       logit_params_for_bad_level_filtering.mat
%           parameters for logit function assigning P_good or P_bad based
%           on the bad level SVM decision boundary distance
%       validation_for_bad_level_filtering.mat
%           validation set on good/bad classification, with overall
%           accuracies used to cap out the logit probabilities returned
%       logit_for_backstep_filter.mat
%           anonymous function for the logit
%       classifiers_for_smart_step_counts.mat
%           set of SVM classifiers with associated logits to make them
%           probabilistic for various step and skip sized from +1 up to +12
%   subroutines:
%       normalizerAC.m
%           function to remove the non-dna-dependent component from the
%           conductance signal, allowing recovery of the continuous curve
%           as a function of DNA position.
%       smartStepCounts.m
%           uses level-to-level continuity information and the associated
%           bad level and step type svm classifiers and logits to calculate
%           at each transition P_step, P_bad, and P_(various skip sizes)
%       levelPredPipeline.m
%           lookup function using a pore model to predict the signal for a
%           given input DNA sequence
%       calibrate_by_iqr.m
%           calibrates two sets of data using a scale and offset to match
%           up their inter quartile ranges
%       shadedErrorBar.m
%           plotting script to generate images with line at mean, shaded
%           region showing error
% 
% 
% FIGURE 3
% ----------------------------------------------------------------------------
% Figure3.m
%   related files:
%       variable_voltage_sequencing_results.mat
%           contains all data quantifying performance of variable voltage
%           sequencing over the reads in this experiment
%       constant_voltage_sequencing_results.mat
%           contains all data quantifying performance of constant voltage
%           sequencing over the reads in this experiment
%       variableVoltageReads/read*section*.mat
%           individual read files for the variable voltage sequencing
%           results. some reads are chimeric--multiple fragments with
%           complementary cutsites re-ligated to create reads where the
%           ground truth is discontiguous portions of the reference genome.
%           these have been partitioned by finding the cutsites in the
%           sequence, and re-aligning the call in between each cut site
%           individually. it is these final individual cut segments that
%           are the "reads" contributing to the final results.
%   subroutines:
%       harvestVariableVoltageSequencingResults.m
%           script to harvest all sequencing results from the variable
%           voltage reads and create the struct used to plot the results
% 
%% FOR CONSTANT VOLTAGE SEQUENCING
% ----------------------------------------------------------------------------
% ----------------------------------------------------------------------------
%
% WRAPPER FOR COMPLETE CONSTANT VOLTAGE PROCESSING
% ----------------------------------------------------------------------------
% processCV.m
%   purpose:
%       wrapper function conducting all processing steps on a constant
%       voltage read: level finding, bad level removal, level
%       recombination, calibration, and sequencing.
%   inputs:
%       read:
%           a constant-voltage read structure, as found in
%           /constantVoltageReads/read_*_section_*.mat
%   outputs:
%       processedread:
%           a read struct like the input read, but now reprocessed
%       sequence:
%           the called sequence for the read
%   subroutines:
%       findLevelsCV.m
%           see below, finds levels in constant voltage data
%       removalFilterCV.m
%           see below, filters out bad levels in constant voltage data
%       recombinationFilterCV.m
%           see below, recombines repeated levels in constant voltage data
%       sequenceCV.m
%           see below, wrapper to sequence constant voltage data
%      
%
% CHANGE POINT DETECTION
% ----------------------------------------------------------------------------
% findLevelsCV.m
%   purpose:
%       finds enzyme steps in constant voltage data
%   inputs:
%       data:
%           current data for a constant voltage read
%   varargin:
%       sensitivity:
%           multiplier for the factor correcting for multiple testing
%           problem. higher values will find fewer levels, smaller values
%           will find more levels
%       minlevellength:
%           minimum length of an allowed level (in data points)
%   outputs:
%       transitions:
%           indices of the transition points
%       features:
%           a 2 x (num_levels) array of [medians ; standard deviations] for
%           each level
%       errors:
%           a 2 x (num_levels) array of the uncertainties in the features
%       stiffnesses:
%           a 1 x (num_levels) cell array of the fisher information
%           matrices for each level
%   related files:
%       cpic_fits_CV.mat
%           contains fits correcting for multiple-testing fallacy in change
%           point detection
%
%
% removalFilterCV.m
%   purpose:
%       finds and removes bad levels from constant voltage data (either
%       gates with low conductance, spikes with high conductance, or very
%       short lived states)
%   inputs:
%       read:
%           read struct for constant voltage read
%   outputs:
%       filteredread:
%           read struct with the bad level removal done
%
%
% RECOMBINATION
% ----------------------------------------------------------------------------
% recombinationFilterCV.m
%   purpose:
%       finds and recombines repeat observations of the same DNA sequence
%       state in constant voltage data. operates similar to the variable
%       voltage recombination filter, using self alignment with a self
%       alignment penalty.
%   inputs:
%       read:
%           a constant voltage read struct, as found in
%           /constantVoltageReads/read_*_section_*.mat
%   varargin:
%       psa:
%           self alignment penalty
%       lookback:
%           how long of a lookback to allow in searching for repeat levels
%       stepprobabilities:
%           probabilities of different enzyme steps (p_step, p_skip,
%           p_skipextend, p_back, p_backextend, p_hold)
%       errorfloor:
%           minimum error to account on any given observation
%   outputs:
%       filteredread:
%           read struct, now with the recombination filter information
%           included
%
%
% CONSTANT VOLTAGE SEQUENCING WRAPPER
% ----------------------------------------------------------------------------
% sequenceCV.m
%   purpose:
%       wrapper to conduct sequencing of constant voltage data. extracts
%       the relevant featuers and stiffnesses, calibrates, loads correct
%       map and transition_info, and calls the sequencer
%   inputs:
%       read:
%           read struct for constant voltage read, as found in
%           /constantVoltageReads/read_*_section_*.mat
%   outputs:
%       calledseq:
%           sequence called for the input read
%   subroutines:
%       calculateSequenceVV.m
%           uses an HMM to determine the maximum likelihood sequence generating
%           an observed set of constant or variable voltage levels.
%   realted files:
%       pore_model_6mer_constant_voltage.mat
%           6-mer pore model for constant voltage with the half-stepping
%           hel308
%       transition_info_hel308_6mer.mat
%           information detailing which kmer model states can transition to
%           which other kmer model states, and with which associated step
%           penalties, assuming a 6-mer pore model and hel308 half stepping
%
% 
%       
%   
% 
%% MISCELLANEOUS
% ----------------------------------------------------------------------------
% ----------------------------------------------------------------------------
% pet28a_sequence.mat
%   a struct containing the sense and antisense strands of the pet28a
%   vector used for sequencing verification experiments, as well as the
%   "complete" sense/antisense sequence, complete = [sense, antisense]
% antisense.m
%   takes as input a 5'->3' DNA sequence and returns the 5'->3'
%   antisense (complement) strand for that sequence.
% alignLetterSequences.m
%   purpose:
%       conducts a simple Needleman-Wunsch style alignment of two input
%       sequences
%   inputs:
%       refseq
%           reference sequence we are aligning to
%       calledseq
%           the sequence we are aligning
%       isglobal
%           boolean of whether to do global (true) or local (false)
%           alignment
%   outputs:
%       aligned
%           3 x length(alignment) character array. 1st row is the reference
%           sequence, 2nd row is the aligned called sequecne, 3rd row
%           represents whether we have a match (+) mismatch (x) insertion
%           (.) or deletion (-)
%       accuracy
%           overall percent accuracy of the alignment, given by N(+) /
%           length(aligned)
% interquartilerange.m
%   purpose: 
%       calculates the interquartile range of a given data set
%   inputs:
%       x
%           data to take the range of
%   outputs:
%       iqr:
%           interquartile range of the data
% jlog.m
%   purpose:
%       computes the jacobian logarithm of a series using an approximation
%       to avoid overflow or roundoff error
%   inputs:
%       d
%           data to take jacobian log of
%   outputs:
%       D
%           resulting jacobian log
%
%   