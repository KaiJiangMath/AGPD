%% sort; first tau, then gamma
function parameters_tot = mysort(parameters_tot)

	%% input:
	%% parameters_tot: all model parameters; size = [nr, 2]
	%% output:
	%% parameters_tot: all model parameters after operating the ascent sort

	%% check and return
	if isempty(parameters_tot)
		fprintf('the input variable is empty in function mysort!\n');
		return;
	end

	%% sort by tau
	[val, ind] = sort(parameters_tot(:,1));
	sort_temp = parameters_tot(ind, :);
	%% the following code lost the data of last line
	sort_temp = [sort_temp; 1000000*ones(1,size(sort_temp,2))];

	%% sort by gamma
	parameters_tot = [];
	ind1 = 1;
	for ij = 1:1:size(sort_temp,1)-1
		if ( abs(sort_temp(ij,1) - sort_temp(ij+1,1)) > 1.0e-8 )
			[val, ind] = sort(sort_temp(ind1:1:ij, 2));
			parameters_tot = [parameters_tot; sort_temp(ind+ind1-1,:)];
			ind1 = ij + 1;
		end
	end
end
