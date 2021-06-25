%% check phase boundary
function boundary_check(model, PATTERN, dim, nr, k2_nc, k2_nr, parameters, check_flag, doc,...
	ini_doc_check, bound_doc_check, calculate, refine_num)

	%% initialization
	[sr, sc] = size(check_flag);

	%% count = [outer_count, bound_count, inner_count]
	%%		pass_count: the number of outer data
	%%		ini_count: the number of adopting new boundary data
	%%		lack_count: the number of inner data
	count = zeros(1, 3);
	
	%% set the range of model parameters
	[taur, gammar] = paraSet(model, calculate);

	%% load energy data about phase diagram
	fstr_pattern = sprintf('%s%s_diagram_hamilton.txt', doc, PATTERN);
	if exist(fstr_pattern) == 2
		ham_matrix = load(fstr_pattern);
	else
		fprintf('\t lack %s\n', fstr_pattern);
		return;
	end
	coordinates = [];

	%% recurrence
	for ij = 1:1:size(parameters,1)
		tau = ignoreNegativeZero(parameters(ij,1));
		gamma = ignoreNegativeZero(parameters(ij,2));

		%% the file path
		fstr1 = sprintf('%s/%s_tau%.6f_gamma%.6f', doc, PATTERN, tau, gamma);
		fstr1_check = sprintf('%s/%s_tau%.6f_gamma%.6f', ini_doc_check, PATTERN, tau, gamma);
		fprintf('\n tau = %f \t gamma = %f\n', tau, gamma);

		%% find [tau, gamma] from ham_matrix
		ind = find_ind(tau, gamma, ham_matrix);
		if ( ind > 0 ) && ( round(ham_matrix(ind,4)) == 0 )	%% symmflag = 0
			count(1) = count(1) + 1;
			fprintf('\t %f\t%f\t outer point\n', tau, gamma);
		elseif ( ind > 0 ) && ( round(ham_matrix(ind,4)) == 2 )	%% lack data
			fprintf('\t %f\t%f\t lack data\n', tau, gamma);
		elseif ( ind > 0 ) && ( round(ham_matrix(ind,4)) == 1 )	%% symmflag = 1
			if exist([fstr1, '.mat']) == 2
				symmflag = 1;
				load([fstr1, '.mat'], 'uCplx', 'rcpBox');
				for i = 1:1:sr	%% all surrounding points
					check_tau = tau + 0.5^refine_num*check_flag(i,1);
					check_gamma = gamma + 0.5^refine_num*check_flag(i,2);
					%% find the positive of tau and gamma in ham_matrix
					ind0 = find_ind(check_tau, check_gamma, ham_matrix);
					if ( ind0 > 0 ) && ( round(ham_matrix(ind0,4)) == 0 ) &&...% symmflag == 0
						( check_tau>=taur(1) && check_tau<=taur(2) ) &&... % the range of tau
						( check_gamma>=gammar(1) && check_gamma<=gammar(2) ) % the range of gamma
						symmflag = 0;
						count(2) = count(2) + 1;
						bound_tau = tau + 0.5^(refine_num+1)*check_flag(i,1); % step/2, space grid refinement
						bound_gamma = gamma + 0.5^(refine_num+1)*check_flag(i,2); % step/2, space grid refinement
						%% find [tau, gamma] from coordinates (should be in a specific range)
						if ( find_ind(bound_tau, bound_gamma, coordinates) == 0 ) &&...
							( bound_tau>=taur(1) && bound_tau<=taur(2) ) &&...
							( bound_gamma>=gammar(1) && bound_gamma<=gammar(2) )
							%% [bound_tau, bound_gamma]: the model parameters waiting computing
							%% [tau, gamma]: the data of the model parameters is used to be the initial value
							coordinates(end+1,:) = [bound_tau, bound_gamma, tau, gamma];
						end
					end
				end
				if ( symmflag == 1 )
					count(3) = count(3) + 1;
					fprintf('\t %f\t%f\t inner point\n', tau, gamma);
				end
			else
				fprintf('\t lack data\n');
			end
		else
			fprintf('\t not in %s\n', fstr_pattern);
		end
	end
	fprintf('\n (The number of) total values: %d\n', sum(count));
	fprintf('\t pass values: %d \t ', count(1));
	fprintf('adopt new initial values: %d \t ', count(2));
	fprintf('lack near values: %d\n', count(3));

	%% sort data; first tau, then gamma
	coordinates = mysort(coordinates);

	%% write parameters data into a specific file
	fid_boundary = sprintf('%sboundary%d.txt', doc, refine_num);
	fidin = fopen(fid_boundary, 'w');
	for ij = 1:1:size(coordinates,1)
		fprintf(fidin, '% .6f\t% .6f\t% .6f\t% .6f\n', coordinates(ij,1), coordinates(ij,2),...
			coordinates(ij,3), coordinates(ij,4));
	end
	fclose(fidin);
end
