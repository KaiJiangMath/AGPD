%%%%%%%%%%%  data:  2021-06-25
function pfc(model, PATTERN, refine_num, i, tau_num, j, gamma_num, calculate, cal_flag)

	%% the code is started
	%clear all; clc; close all
	set(0, 'DefaultFigureVisible', 'off')
	%calculate = 'reflag', 'data_check', 'diagram_check' or others;
	%cal_flag = 'yes' or 'no';
	fig_flag = 'yes'; % plotting figures; 'yes' or 'no'

%	model = 'lp';
%	PATTERN = 'LQ6';
%	refine_num = '0';
%	i = '0';
%	j = '0';
%	tau_num = '1';
%	gamma_num = '1';
%	calculate = 'compute';
%	cal_flag = 'yes';

	%profile on 
	refine_num = round(str2num(refine_num));
	scaleQ = 1.0;

	%% folder for each pattern
	doc = sprintf('%s_results/%s/', model, PATTERN);
	if exist(doc) == 0
		mkdir(doc);
	end

	%% add the function path
	addpath(genpath('src/'));

	%% load the initial parameters
	[taur_tot, gammar_tot] = paraSet(model, calculate);

	%% split region
	i = round(str2num(i));		tau_num = round(str2num(tau_num));
	j = round(str2num(j));		gamma_num = round(str2num(gamma_num));
	fprintf('tau split: \t')
	taur = split_region(taur_tot, i, tau_num)
	fprintf('gamma split: \t')
	gammar = split_region(gammar_tot, j, gamma_num)
	%% the step of tau and gamma
	dtau = taur(3);
	dg = gammar(3);

	%% check_flag for searching surrounding points
	check_flag = [	-1,	0; 
					1,	0; 
					0,	-1; 
					0,	1;
					-1,	-1;
					-1,	1;
					1, -1;
					1,	1];
	check_flag = [dtau*check_flag(:,1), dg*check_flag(:,2)];

	%% set the model parameters
	%% model parameters
	parameters = [];
	for tauj = taur(1):dtau:taur(2)
		for gammaj = gammar(1):dg:gammar(2)
			parameters(end+1,:) = [tauj, gammaj, 1000000, 1000000];
		end
	end
	%% add boundary model parameters
	for rj = 0:1:refine_num
		fid_boundary = sprintf('%sboundary%d.txt', doc, rj);
		if ( exist(fid_boundary) == 2 )
			temp = load(fid_boundary);
			%% split the model parameters to accelerate computations
			index = split_region([1, size(temp,1), 1], i*gamma_num+j, tau_num*gamma_num);
			temp = temp(index(1):1:index(2), :);
			%% all model parameters
			parameters = [parameters; temp];
		end
	end
	%% ignore negative zero
	for ij = 1:1:size(parameters,1)
		for pj = 1:1:size(parameters,2)
			parameters(ij,pj) = ignoreNegativeZero(parameters(ij,pj));
		end
	end


	%% spatial discretization
	grid = 'coarse';
	%grid = 'accurate';
	if strcmp(grid, 'coarse')
		if strcmp(PATTERN, 'lam') || strcmp(PATTERN, 'hex')
			nc = 24*[1,1];		nr = 48*[1,1];
		elseif strcmp(PATTERN, 'sigma')	
			nc = 64*[2,2,1];	nr = 64*[2,2,1];
		elseif strcmp(PATTERN, '10fold') || strcmp(PATTERN, '12fold')
			nc = 20*[1,1,1,1];	nr = 20*[1,1,1,1];
		elseif strcmp(PATTERN, 'LQ6') || strcmp(PATTERN, 'LQS6') || strcmp(PATTERN, 'C3') ||...
		   	strcmp(PATTERN, '12i6o') || strcmp(PATTERN, '8i10o') || strcmp(PATTERN, 'sq') ||...
			strcmp(PATTERN, 'squ') || strcmp(PATTERN, 'sqv') || strcmp(PATTERN, 'sqw') ||...
			strcmp(PATTERN, 'Ls')
			nc = 20*[1,1,1,1];	nr = 20*[1,1,1,1];
		else
			nc = 24*[1,1,1];	nr = 48*[1,1,1];
		end
	elseif strcmp(grid, 'accurate')
		if strcmp(PATTERN, 'lam') || strcmp(PATTERN, 'hex')
			nc = 48*[1,1];		nr = 64*[1,1];
		elseif strcmp(PATTERN, 'sigma')	
			nc = 128*[2,2,1];	nr = 128*[2,2,1];
		elseif strcmp(PATTERN, '10fold') || strcmp(PATTERN, '12fold')
			nc = 24*[1,1,1,1];	nr = 24*[1,1,1,1];
		elseif strcmp(PATTERN, 'LQ6') || strcmp(PATTERN, 'LQS6') || strcmp(PATTERN, 'C3') ||...
		   	strcmp(PATTERN, '12i6o') || strcmp(PATTERN, '8i10o') || strcmp(PATTERN, 'sq') ||...
			strcmp(PATTERN, 'squ') || strcmp(PATTERN, 'sqv') || strcmp(PATTERN, 'sqw') ||...
			strcmp(PATTERN, 'Ls')
			nc = 24*[1,1,1,1];	nr = 24*[1,1,1,1];
		else
			nc = 48*[1,1,1];	nr = 64*[1,1,1];
		end
	else
		fprintf('WARNING: without setting grid!\n');
	end

	%% check file just for copying figures;
	if strcmp(calculate, 'data_check') || strcmp(calculate, 'diagram_check') ||...
		strcmp(calculate, 'boundary')
		doc_check = sprintf('%scheck/', doc);
		if exist(doc_check) == 0
			mkdir(doc_check);
		end
		ini_doc_check = sprintf('%sini_check/', doc);
		if ~strcmp(cal_flag, 'yes') && exist(ini_doc_check) == 0
			mkdir(ini_doc_check);
		end
	end


	%% initial value and box;
	[uinit, rcpBox, dim] = IniCplxBox(PATTERN, nc, scaleQ);

	tic 
	projmat = getprojmat(PATTERN, dim(1), dim(2));
	k2_nc = obtGsquare(nc, rcpBox, dim, projmat);
	if max(abs(nc-nr)) < 1.0e-8
		k2_nr = k2_nc;
	else
		k2_nr = obtGsquare(nr, rcpBox, dim, projmat);
	end
	toc

	%% plot initial value;
	dirBox = getDualBox(rcpBox, scaleQ)
	rcpBox
	if ( strcmp(fig_flag, 'yes') )
		init_fig = sprintf('%sinit', doc);
		if exist([init_fig, '.png']) == 0
			plotPhase(PATTERN, uinit, dirBox, dim, rcpBox, init_fig);
			save([init_fig, '.mat'], 'uinit', 'dirBox', 'dim', 'rcpBox');
		end
		if ( strcmp(calculate, 'data_check') || strcmp(calculate, 'diagram_check') ||...
			strcmp(calculate, 'boundary') ) && exist([doc_check, 'init.png']) == 0
			copyfile([init_fig, '.png'], [doc_check, 'init.png']);
		end
	end


	fprintf('\n\n\n=================> %s <===============\n\n\n', PATTERN)
	fname = sprintf('%s%s_hamilton.txt', doc, PATTERN);
	%% check based on the existing data;
	if strcmp(calculate, 'reflag')
		if exist(fname) == 0
			fprintf('WARNING: check %s please!\n', fname);
			return;
		else
			ham_matrix = load(fname);
			for ij = 1:1:size(parameters,1)
				tau = parameters(ij,1);
				gamma = parameters(ij,2);
				fstr1 = sprintf('%s/%s_tau%.6f_gamma%.6f', doc, PATTERN, tau, gamma);
				fprintf('\n tau = %f \t gamma = %f\n', tau, gamma);
				%%	
				ind = find_ind(tau, gamma, ham_matrix);
				if ( ind > 0 )
					hamilton = ham_matrix(ind, end-1);
					symmflag = ham_matrix(ind, end);
					if ( symmflag ~= 2 )
						if ( hamilton > -5.0e-8 && symmflag == 1)
							symmflag = 0;
						end
						save([fstr1, '.mat'], 'symmflag', '-append');
					end
				else
					fprintf('\t%f\t%f\t lack in %s\n', tau, gamma, fname);
				end
			end
		end
	elseif strcmp(calculate, 'data_check')
		data_check(model, PATTERN, dim, nr, k2_nc, k2_nr, parameters, check_flag, doc,...
			doc_check, ini_doc_check, calculate, cal_flag, refine_num, fig_flag);
	elseif strcmp(calculate, 'diagram_check')
		diagram_check(model, PATTERN, dim, nr, k2_nc, k2_nr, parameters, check_flag, doc,...
			doc_check, ini_doc_check, calculate, cal_flag, refine_num, fig_flag);
	elseif strcmp(calculate, 'boundary')
		boundary_check(model, PATTERN, dim, nr, k2_nc, k2_nr, parameters, check_flag, doc,...
			ini_doc_check, doc_check, calculate, refine_num);
	elseif strcmp(calculate, 'compute')
		compute(model, PATTERN, dim, uinit, nr, k2_nc, k2_nr, parameters, doc,...
			rcpBox, calculate, cal_flag, fig_flag);
	elseif strcmp(calculate, 'tidy')
		tidy(fname, model, PATTERN, parameters, doc, refine_num, check_flag, cal_flag);
	else
		fprintf('check %s please\n', calculate);
	end

	fprintf('%s \t %s \t %s \t %s\n', model, PATTERN, calculate, cal_flag);
	fprintf('the refine number: %d\n', refine_num);
	fprintf('tau: start: %.6f \t end: %.6f \t step: %.6f \t split: %d\n',...
		taur(1), taur(2), taur(3), tau_num);
	fprintf('gamma: start: %.6f \t end: %.6f \t step: %.6f \t split: %d\n',...
		gammar(1), gammar(2), gammar(3), gamma_num);

	%% remove the function path
	rmpath(genpath('src/'));

	%profile viewer

	%% the code is finished
	finish_file = sprintf('finish/%s-[%d-%d]-[%d-%d].txt', PATTERN, i, tau_num, j, gamma_num);
	finish = fopen(finish_file, 'w');
	fprintf(finish, 'finish');
	fclose(finish);
end
