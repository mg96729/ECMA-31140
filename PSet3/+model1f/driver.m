%
% Status : main Dynare file
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

if isoctave || matlab_ver_less_than('8.6')
    clear all
else
    clearvars -global
    clear_persistent_variables(fileparts(which('dynare')), false)
end
tic0 = tic;
% Define global variables.
global M_ options_ oo_ estim_params_ bayestopt_ dataset_ dataset_info estimation_info ys0_ ex0_
options_ = [];
M_.fname = 'model1f';
M_.dynare_version = '5.5';
oo_.dynare_version = '5.5';
options_.dynare_version = '5.5';
%
% Some global variables initialization
%
global_initialization;
M_.exo_names = {};
M_.exo_names_tex = {};
M_.exo_names_long = {};
M_.endo_names = cell(2,1);
M_.endo_names_tex = cell(2,1);
M_.endo_names_long = cell(2,1);
M_.endo_names(1) = {'k'};
M_.endo_names_tex(1) = {'k'};
M_.endo_names_long(1) = {'k'};
M_.endo_names(2) = {'c'};
M_.endo_names_tex(2) = {'c'};
M_.endo_names_long(2) = {'c'};
M_.endo_partitions = struct();
M_.param_names = cell(7,1);
M_.param_names_tex = cell(7,1);
M_.param_names_long = cell(7,1);
M_.param_names(1) = {'alpha'};
M_.param_names_tex(1) = {'alpha'};
M_.param_names_long(1) = {'alpha'};
M_.param_names(2) = {'beta'};
M_.param_names_tex(2) = {'beta'};
M_.param_names_long(2) = {'beta'};
M_.param_names(3) = {'delta'};
M_.param_names_tex(3) = {'delta'};
M_.param_names_long(3) = {'delta'};
M_.param_names(4) = {'sigma'};
M_.param_names_tex(4) = {'sigma'};
M_.param_names_long(4) = {'sigma'};
M_.param_names(5) = {'k_s'};
M_.param_names_tex(5) = {'k\_s'};
M_.param_names_long(5) = {'k_s'};
M_.param_names(6) = {'c_s'};
M_.param_names_tex(6) = {'c\_s'};
M_.param_names_long(6) = {'c_s'};
M_.param_names(7) = {'k_0'};
M_.param_names_tex(7) = {'k\_0'};
M_.param_names_long(7) = {'k_0'};
M_.param_partitions = struct();
M_.exo_det_nbr = 0;
M_.exo_nbr = 0;
M_.endo_nbr = 2;
M_.param_nbr = 7;
M_.orig_endo_nbr = 2;
M_.aux_vars = [];
M_ = setup_solvers(M_);
M_.Sigma_e = zeros(0, 0);
M_.Correlation_matrix = eye(0, 0);
M_.H = 0;
M_.Correlation_matrix_ME = 1;
M_.sigma_e_is_diagonal = true;
M_.det_shocks = [];
M_.surprise_shocks = [];
M_.heteroskedastic_shocks.Qvalue_orig = [];
M_.heteroskedastic_shocks.Qscale_orig = [];
options_.linear = false;
options_.block = false;
options_.bytecode = false;
options_.use_dll = false;
M_.orig_eq_nbr = 2;
M_.eq_nbr = 2;
M_.ramsey_eq_nbr = 0;
M_.set_auxiliary_variables = exist(['./+' M_.fname '/set_auxiliary_variables.m'], 'file') == 2;
M_.epilogue_names = {};
M_.epilogue_var_list_ = {};
M_.orig_maximum_endo_lag = 1;
M_.orig_maximum_endo_lead = 1;
M_.orig_maximum_exo_lag = 0;
M_.orig_maximum_exo_lead = 0;
M_.orig_maximum_exo_det_lag = 0;
M_.orig_maximum_exo_det_lead = 0;
M_.orig_maximum_lag = 1;
M_.orig_maximum_lead = 1;
M_.orig_maximum_lag_with_diffs_expanded = 1;
M_.lead_lag_incidence = [
 1 2 0;
 0 3 4;]';
M_.nstatic = 0;
M_.nfwrd   = 1;
M_.npred   = 1;
M_.nboth   = 0;
M_.nsfwrd   = 1;
M_.nspred   = 1;
M_.ndynamic   = 2;
M_.dynamic_tmp_nbr = [2; 0; 0; 0; ];
M_.model_local_variables_dynamic_tt_idxs = {
};
M_.equations_tags = {
  1 , 'name' , '1' ;
  2 , 'name' , '2' ;
};
M_.mapping.k.eqidx = [1 2 ];
M_.mapping.c.eqidx = [1 2 ];
M_.static_and_dynamic_models_differ = false;
M_.has_external_function = false;
M_.state_var = [1 ];
M_.exo_names_orig_ord = [1:0];
M_.maximum_lag = 1;
M_.maximum_lead = 1;
M_.maximum_endo_lag = 1;
M_.maximum_endo_lead = 1;
oo_.steady_state = zeros(2, 1);
M_.maximum_exo_lag = 0;
M_.maximum_exo_lead = 0;
oo_.exo_steady_state = zeros(0, 1);
M_.params = NaN(7, 1);
M_.endo_trends = struct('deflator', cell(2, 1), 'log_deflator', cell(2, 1), 'growth_factor', cell(2, 1), 'log_growth_factor', cell(2, 1));
M_.NNZDerivatives = [6; -1; -1; ];
M_.static_tmp_nbr = [2; 1; 0; 0; ];
M_.model_local_variables_static_tt_idxs = {
};
M_.params(1) = 0.36;
alpha = M_.params(1);
M_.params(2) = 0.95;
beta = M_.params(2);
M_.params(3) = 0.025;
delta = M_.params(3);
M_.params(4) = 2;
sigma = M_.params(4);
M_.params(5) = (M_.params(1)*M_.params(2)/(1+M_.params(2)*(M_.params(3)-1)))^(1/(1-M_.params(1)));
k_s = M_.params(5);
M_.params(6) = M_.params(5)^M_.params(1)-M_.params(3)*M_.params(5);
c_s = M_.params(6);
M_.params(7) = M_.params(5)*0.1;
k_0 = M_.params(7);
%
% INITVAL instructions
%
options_.initval_file = false;
oo_.steady_state(1) = M_.params(7);
if M_.exo_nbr > 0
	oo_.exo_simul = ones(M_.maximum_lag,1)*oo_.exo_steady_state';
end
if M_.exo_det_nbr > 0
	oo_.exo_det_simul = ones(M_.maximum_lag,1)*oo_.exo_det_steady_state';
end
%
% ENDVAL instructions
%
ys0_= oo_.steady_state;
ex0_ = oo_.exo_steady_state;
oo_.steady_state(1) = M_.params(5);
oo_.steady_state(2) = M_.params(6);
resid;
options_.periods = 200;
perfect_foresight_setup;
perfect_foresight_solver;
var_list_ = {'k'};
rplot(var_list_);


oo_.time = toc(tic0);
disp(['Total computing time : ' dynsec2hms(oo_.time) ]);
if ~exist([M_.dname filesep 'Output'],'dir')
    mkdir(M_.dname,'Output');
end
save([M_.dname filesep 'Output' filesep 'model1f_results.mat'], 'oo_', 'M_', 'options_');
if exist('estim_params_', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'model1f_results.mat'], 'estim_params_', '-append');
end
if exist('bayestopt_', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'model1f_results.mat'], 'bayestopt_', '-append');
end
if exist('dataset_', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'model1f_results.mat'], 'dataset_', '-append');
end
if exist('estimation_info', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'model1f_results.mat'], 'estimation_info', '-append');
end
if exist('dataset_info', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'model1f_results.mat'], 'dataset_info', '-append');
end
if exist('oo_recursive_', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'model1f_results.mat'], 'oo_recursive_', '-append');
end
if ~isempty(lastwarn)
  disp('Note: warning(s) encountered in MATLAB/Octave code')
end
