%%%%%%%%%%%  Computing ordered  on phase field model by CH equation
%%%%%%%%%%%  Projection method with P=I is used in this code
%%%%%%%%%%%  The semi-implicit scheme is the primary algorithm to minimize the problem
%%%%%%%%%%%  At the same time, the Nesterov method is used to accelerate the optimization algorithm
%%%%%%%%%%%  The time step size is estimated by mean eigevalue of Hessian matrix
%%%%%%%%%%%  data:  2020-03-03

function [hamilton, symmflag] = pfcmodel(model, PATTERN, dim, uCplx, nr, k2_nc,...
	k2_nr, parameter, rcpBox, calculate, fname, fig_flag)
	
    format long;
	close all;

	if max(abs(size(k2_nc) - size(k2_nr))) == 0
		k2flag = 0;
	else
		k2flag = 1;
	end

	ncpt = size(uCplx);

	uCplx_cs = uCplx;
	if dim(1) == 1 uCplx_rf = zeros(1, nr);
	else uCplx_rf = zeros(nr);
	end

	%% the type of flow: AC or CH
	flow_type = 'AC';
%	flow_type = 'CH';

	if strcmp(PATTERN, 'hex')
		type = 'hex';
	elseif strcmp(PATTERN, 'sigma')
		type = 'rectangular';
	else
		type = 'cube';
	end

	scaleQ = 1;
	DimCpt = dim(1);
	DimPhy = dim(2);

	%%%%%% model systems %%%%%%%
	sysPmts(1) = 1.0;   %%%  scale 
	if strcmp(model, 'lp')
		sysPmts(2) = 2*cos(pi/12);   %%%  another scale 
	else
		sysPmts(2) = 0;   %%%  another scale 
	end
	sysPmts(3) = 1.0;   %%%  interaction term
%	sysPmts(3) = 0.65^2;   %%%  interaction term
	if strcmp(model, 'lp')
		sysPmts(4) = parameter(1);	%%%  2nd term
		sysPmts(5) = 2*parameter(2); %%%  3rd term % without factorial
		sysPmts(6) = 6;   %%%  4th term % without factorial
	else
		sysPmts(4) = parameter(1);   %%%  2nd term
		sysPmts(5) = parameter(2);   %%%  3rd term
		sysPmts(6) = 1;   %%%  4th term
	end

%	tstep = 1.0;
	tstep = 0.1;
	tmin = 1.0e-2;
	tmax = 1.0;

	t0 = tstep;
	diffham = inf;
	hamOld = inf;
	hamilton = inf;
	iterator = 0;
	err = 1.0;
	itMax = 5000;
	tol_cs = 5.0e-4;
	tol_rf = 1.0e-6;
	if ( k2flag == 0 ); tol_cs = tol_rf; end

	cmpPmts(1) = tol_cs;
	cmpPmts(2) = itMax;
	cmpPmts(3) = tstep;
	cmpPmts(4) = tmin;
	cmpPmts(5) = tmax;
	
	EPS = 1.0;
	cout = 0;
	if strcmp(model, 'lp')
		cout_max = 0;
	else
		cout_max = 20;
	end

	tic
	%%% fixed box
	if strcmp(calculate, 'data_check') || strcmp(calculate, 'diagram_check') ||...
		strcmp(calculate, 'boundary')
		uCplx_rf = uCplx;
	elseif ( k2flag == 0 )
		uCplx_rf = uCplx;
	else
		[ham0, uCplx_cs] = pfcapg(model, flow_type, PATTERN,...
			k2_nc, sysPmts, dim, cmpPmts, rcpBox, uCplx, fname);
		uCplx_rf = mymap(uCplx_cs, uCplx_rf, 'cs2rf');
		fprintf('coarse grid to refine grid\n');
	end
	cmpPmts(1) = tol_rf;
	[ham0, uCplx_rf] = pfcapg(model, flow_type, PATTERN,...
		k2_nr, sysPmts, dim, cmpPmts, rcpBox, uCplx_rf, fname);

	while (EPS > 1.0e-6 && cout < cout_max)
		cout = cout+1;

		%% show box in directly space and reciprocal space
		[rcpBox, k2] = pfcbox(model, PATTERN, sysPmts, dim, uCplx_rf, rcpBox, 1.0e-5, type);
		fprintf('\n\n');
		[ham1, uCplx_rf] = pfcapg(model, flow_type, PATTERN,...
			k2, sysPmts, dim, cmpPmts, rcpBox, uCplx_rf, fname);
		EPS = abs(ham1 - ham0);
		ham0 = ham1;
		fprintf('cout %d : hamilton = %.15e\t EPS = %.15e\n\n\n', cout, ham0, EPS);

	end
	toc

	%% plot
	dirBox = getDualBox(rcpBox, scaleQ)
	rcpBox
	if strcmp(fig_flag, 'yes')
		plotPhase(PATTERN, uCplx_rf, dirBox, dim, rcpBox, fname);
	end

	%% identify
	symmflag = phase_identify(PATTERN, uCplx_rf, nr, dim, fig_flag, fname);

	hamilton = ham0
	ncpt = nr;		uCplx = uCplx_rf;
    save([fname, '.mat'], 'uCplx', 'ncpt', 'rcpBox', 'dim', 'hamilton', 'symmflag');
end


