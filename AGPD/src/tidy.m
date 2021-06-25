%% clear up all data
function tidy(fname, model, PATTERN, parameters, doc, refine_num, check_flag, cal_flag)

	%% input:
	%% cal_flag: a special parameter for the function tidy
	%%			here we give two choices: 'parameter' and 'hamilton'

	%% set the range of model parameters
	[taur, gammar] = paraSet(model, 'tidy');
	[sr, sc] = size(check_flag);

	%% different case
	if strcmp(cal_flag, 'parameters')
		%% clear up the waiting computing parameters of various patterns
		fid_boundary = sprintf('%sboundary%d.txt', doc, refine_num);
		parameters_tot = [];
		if ( exist(fid_boundary) == 2 )
			parameters = load(fid_boundary); % load existing data
			for ij = 1:1:size(parameters,1)
				%% if the point is not in parameters_tot
				if ( find_ind(parameters(ij,1), parameters(ij,2), parameters_tot) == 0 )
					parameters_tot(end+1,:) = parameters(ij,:); % all parameters
					for i1 = 1:1:sr	%% all surrounding points
						check_tau = parameters(ij,1) + 0.5^(refine_num+1)*check_flag(i1,1); % step/2
						check_gamma = parameters(ij,2) + 0.5^(refine_num+1)*check_flag(i1,2); % step/2
						if ( find_ind(check_tau, check_gamma, parameters_tot) == 0 ) &&...
							( check_tau>=taur(1) && check_tau<=taur(2) ) &&... % the range of tau
							( check_gamma>=gammar(1) && check_gamma<=gammar(2) ) % the range of gamma
							parameters_tot(end+1,:) = [check_tau, check_gamma, parameters(ij,1:1:2)]; % add extra points
						end
					end
				end
			end
		else
			fprintf('WARNING: lack %s\n', fid_boundary);
		end
		if isempty(parameters_tot)
			return;
		end
		%% sort; first tau, then gamma
		parameters_tot = mysort(parameters_tot);
		%% clear up all waiting computing parameters
		fid_boundary = sprintf('%sboundary%d.txt', doc, refine_num);
		fidin = fopen(fid_boundary, 'w');
		for ij = 1:1:size(parameters_tot,1)
			fprintf(fidin, '% .6f\t% .6f\t% .6f\t% .6f\n', parameters_tot(ij,1),...
				parameters_tot(ij,2), parameters_tot(ij,3), parameters_tot(ij,4));
		end
		fclose(fidin);
	elseif strcmp(cal_flag, 'hamilton')
		%% how many parameters are computed
		exist_count = 0;
		eham = fopen(fname, 'w');
		for ij = 1:1:size(parameters,1)
			tau = ignoreNegativeZero(parameters(ij,1));
			gamma = ignoreNegativeZero(parameters(ij,2));
			%% the file path
			fstr1 = sprintf('%s/%s_tau%.6f_gamma%.6f', doc, PATTERN, tau, gamma);
			fprintf('\n tau = %f \t gamma = %f\n', tau, gamma);
			if exist([fstr1, '.mat']) == 0		%% without data
				fprintf(eham, '% .6f\t% .6f\t % .6e\t%d\n', tau, gamma, 0, 2);	%% 2: lack data
				fprintf('\t % .6f\t% .6f\t lack data\n', tau, gamma);
			else
				exist_count = exist_count + 1;
				load([fstr1, '.mat'], 'hamilton', 'symmflag');
				if isnan(hamilton)	hamilton = 0; end
				fprintf(eham, '% .6f\t% .6f\t % .6e\t%d\n', tau, gamma, hamilton, symmflag);
			end
		end
		fclose(eham);
		fprintf('\n (The number of) total values: %d\n', size(parameters,1));
		fprintf('exist values: %d \t noexist values: %d\n', exist_count, size(parameters,1)-exist_count);
	else
		fprintf('WARNING: the parameter "cal_flag" is wrong in function tidy\n');
	end
end
