%% identify structures by aHash algorithms;
function symmflag = phase_identify(PATTERN, uCplx, ncpt, dim, fig_flag, fname);

	%% input:
	%% PATTERN: the computing candidate structures;
	%% uCplx: the order parameter in reciprocal space;
	%% ncpt: the spatial discretization;
	%% dim: the computational dimensionality;
	%% fig_flag: 'yes' or 'no' to determine whether figures are plotted;
	%% fname: the figure saved path;
	%% output:
	%% symmflag: the flag of symmetry;
	%%		1: maintaining symmetry; 0: breaking symmetry;

	%% set the number of identifying samples;
	idnum = 0; % the default value;
	switch PATTERN
	case 'hcp'
		idnum = 1;
	case 'bcc'
		idnum = 5;
	case 'fcc'
		idnum = 4;
	case 'gyroid'
		idnum = 2;
	case 'A15'
		idnum = 3;
	case 'sigma'
		idnum = 1;
	case 'fddd'
		idnum = 1;
	case 'fdddCube'
		idnum = 1;
	end

	%% check structures' symmetry according to the value of "idnum";
	if ( idnum == 0 )
		uc0 = obtainIniVals(PATTERN, ncpt, dim);
		symmflag = detsymmetry(uc0, uCplx);
	else
		if strcmp(fig_flag, 'yes')
			for j = 1:1:idnum
				%% load the identifying sample to check the computing phase;
				file_ref = sprintf('identify/%s/%d.png', PATTERN, j);
				[similarity, symmflag] = direct_identify(PATTERN, [fname, '.png'], file_ref);
				if ( symmflag == 1 ) break; end
			end
		else
			symmflag = 1;
		end
	end

	%% give annotation;
	if ( symmflag == 1 )
        fprintf('--> Maintain Symmetry ...\n');
	else
		fprintf('--> WARNING: Symmetry Breaking ...\n');
	end
end
