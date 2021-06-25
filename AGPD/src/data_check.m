%% check data
%% compute data from a better initial value
function data_check(model, PATTERN, dim, nr, k2_nc, k2_nr, parameters, check_flag, doc,...
	doc_check, ini_doc_check, calculate, cal_flag, refine_num, fig_flag)

	%% initialization
	[sr, sc] = size(check_flag);

	%% count = [pass_count, ini_count, lack_count]
	%%		pass_count: the number of pass data
	%%		ini_count: the number of adopting new initial data
	%%		lack_count: the number of the data which lacks near values
	count = zeros(1, 3);

	for ij = 1:1:size(parameters,1)
		tau = parameters(ij,1);
		gamma = parameters(ij,2);
		%% file path
		fstr1 = sprintf('%s/%s_tau%.6f_gamma%.6f', doc, PATTERN, tau, gamma);
		fstr1_check = sprintf('%s/%s_tau%.6f_gamma%.6f', doc_check, PATTERN, tau, gamma);
		fstr2_check = sprintf('%s/%s_tau%.6f_gamma%.6f', ini_doc_check, PATTERN, tau, gamma);
		fprintf('\n tau = %f \t gamma = %f\n', tau, gamma);
		%%
		if exist([fstr1, '.mat']) == 2
			load([fstr1, '.mat'], 'symmflag');
			if ( symmflag == 1 )
				count(1) = count(1) + 1;
				fprintf('\t %f\t%f\t pass check\n', tau, gamma);
			else
				for i = 1:1:sr
					ini_tau = tau + 0.5^refine_num*check_flag(i,1);
					ini_gamma = gamma + 0.5^refine_num*check_flag(i,2);
					ini_fstr1 = sprintf('%s/%s_tau%.6f_gamma%.6f',...
						doc, PATTERN, ini_tau, ini_gamma);
					ini_fstr1_check = sprintf('%s/%s_tau%.6f_gamma%.6f',...
						ini_doc_check, PATTERN, ini_tau, ini_gamma);
					if ( exist([ini_fstr1, '.mat']) == 2 )
						load([ini_fstr1, '.mat'], 'symmflag');
						if ( symmflag == 1 )
							count(2) = count(2) + 1;
							fprintf('\t %f\t%f\t initial data\n', ini_tau, ini_gamma);
							if strcmp(cal_flag, 'yes')
								load([ini_fstr1, '.mat'], 'uCplx', 'rcpBox');
								%% obtain ksquare for using the existing data
								if abs(size(uCplx)-nr) > 1.0e-8
									nr = size(uCplx);
									projmat = getprojmat(PATTERN, dim(1), dim(2));
									k2_nr = obtGsquare(nr, rcpBox', dim, projmat');
									k2_nc = k2_nr;
								end
								[hamilton, symmflag] = pfcmodel(model, PATTERN, dim, uCplx,...
									nr, k2_nc, k2_nr, parameters(ij,:), rcpBox, calculate, fstr1, fig_flag);
								if strcmp(fig_flag, 'yes')
									copyfile([fstr1, '.png'], [fstr1_check, '.png']);
								end
							else
								copyfile([ini_fstr1, '.mat'], [ini_fstr1_check, '.mat']);
								copyfile([fstr1, '.mat'], [fstr2_check, '.mat']);
								if strcmp(fig_flag, 'yes')
									copyfile([ini_fstr1, '.png'], [ini_fstr1_check, '.png']);
									copyfile([fstr1, '.png'], [fstr2_check, '.png']);
								end
							end
							break;
						end
					end
				end
				if ( symmflag == 0 )
					count(3) = count(3) + 1;
					fprintf('\t %f\t%f\t lack near data\n', tau, gamma);
				end
			end
		else
			fprintf('\t %f\t%f\t lack in hamilton.txt\n', tau, gamma);
		end
	end
	fprintf('\n (The number of) total values: %d\n', sum(count));
	fprintf('\t pass values: %d \t ', count(1));
	fprintf('adopt new initial values: %d \t ', count(2));
	fprintf('lack near values: %d\n', count(3));
end
