%% compute data
function compute(model, PATTERN, dim, uinit, nr, k2_nc, k2_nr, parameters, doc,...
	rcpBox, calculate, cal_flag, fig_flag)

	exist_count = 0;
	for ij = 1:1:size(parameters,1)
		tau = parameters(ij,1);
		gamma = parameters(ij,2);
		%% the file path
		fstr1 = sprintf('%s/%s_tau%.6f_gamma%.6f', doc, PATTERN, tau, gamma);
		fprintf('\n tau = %f \t gamma = %f\n', tau, gamma);
		if exist([fstr1, '.mat']) == 0		%% without data
			if strcmp(cal_flag, 'yes')		%% compute
				init_file = sprintf('%s/%s_tau%.6f_gamma%.6f', doc, PATTERN, parameters(ij,3), parameters(ij,4));
				if exist(init_file) == 2
					load(init_file);
				else
					uCplx = uinit;
				end
				[hamilton, symmflag] = pfcmodel(model, PATTERN, dim, uCplx, nr, k2_nc,...
					k2_nr, parameters(ij,:), rcpBox, calculate, fstr1, fig_flag);
				if isnan(hamilton) hamilton = 0; end
			else
				fprintf('\t %f\t%f\t lack data\n', tau, gamma);
			end
		else
			exist_count = exist_count + 1;
			load([fstr1, '.mat'], 'hamilton', 'symmflag');
			if isnan(hamilton)	hamilton = 0; end
			fprintf('hamilton: %.6e\n', hamilton);
		end
	end
	fprintf('\n (The number of) total values: %d\n', size(parameters,1));
	fprintf('exist values: %d \t noexist values: %d\n', exist_count, size(parameters,1)-exist_count);
end
