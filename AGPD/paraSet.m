%% set the model parameters
function [taur_tot, gammar_tot] = paraSet(model, calculate)

	%% set the model parameters (tau: coefficient before 2nd term, gamma: 3rd term)
	%% the energy penalty factor c is given in the file named pfcmodel.m (rescaling as 1)

	if strcmp(model, 'lb')
		%% space step
		dtau = 0.05;		dg = 0.1;
		if strcmp(calculate, 'tidy') || strcmp(calculate, 'reflag')
			%% there should be complete data for plotting phase diagram
			taur_tot = [-0.4, 0.4, dtau];
			gammar_tot = [0.0, 2.0, dg];
%			taur_tot = [-0.2, 0.3, dtau];
%			gammar_tot = [0.8, 1.4, dg];
		else
			%% data can be split to accelerate for one certain pattern
			taur_tot = [-0.4, 0.4, dtau];
			gammar_tot = [0.0, 2.0, dg];
%			taur_tot = [-0.2, 0.3, dtau];
%			gammar_tot = [0.8, 1.4, dg];
		end
	elseif strcmp(model, 'lp')
		%% there are no factorials in the nonlinear part of the LP free energy functional
		%% just set tau = epsilon (coefficient before 2nd term), gamma = alpha (3rd term)
		%% space step
		dtau = 0.002;		dg = 0.05;
		if strcmp(calculate, 'tidy') || strcmp(calculate, 'reflag')
			%% there should be complete data for plotting phase diagram
			taur_tot = [-0.01, 0.05, dtau];
			gammar_tot = [0.0, 1.0, dg];
%			taur_tot = [-0.01, 0.10, dtau];
%			gammar_tot = [0.0, 2.0, dg];
		else
			%% data can be split to accelerate for one certain pattern
			taur_tot = [-0.01, 0.05, dtau];
			gammar_tot = [0.0, 1.0, dg];
%			taur_tot = [-0.01, 0.10, dtau];
%			gammar_tot = [0.0, 2.0, dg];
		end
	elseif strcmp(model, 'ok')
		%% space step
		dtau = 0.05;		dg = 0.1;
		if strcmp(calculate, 'tidy') || strcmp(calculate, 'reflag')
			%% there should be complete data for plotting phase diagram
			taur_tot = [-0.4, 0.4, dtau];
			gammar_tot = [0.0, 2.0, dg];
		else
			%% data can be split to accelerate for one certain pattern
			taur_tot = [-0.4, 0.4, dtau];
			gammar_tot = [0.0, 2.0, dg];
		end
	elseif strcmp(model, 'leibler')
		%% space step
		dtau = 0.05;		dg = 0.1;
		if strcmp(calculate, 'tidy') || strcmp(calculate, 'reflag')
			%% there should be complete data for plotting phase diagram
			taur_tot = [-0.4, 0.4, dtau];
			gammar_tot = [0.0, 2.0, dg];
		else
			%% data can be split to accelerate for one certain pattern
			taur_tot = [-0.4, 0.4, dtau];
			gammar_tot = [0.0, 2.0, dg];
		end
	else
		%% space step
		dtau = 0.05;		dg = 0.1;
		if strcmp(calculate, 'tidy') || strcmp(calculate, 'reflag')
			%% there should be complete data for plotting phase diagram
			taur_tot = [-0.4, 0.4, dtau];
			gammar_tot = [0.0, 2.0, dg];
		else
			%% data can be split to accelerate for one certain pattern
			taur_tot = [-0.4, 0.4, dtau];
			gammar_tot = [0.0, 2.0, dg];
		end
	end

	%% test code
%	dtau = 0.01;		dg = 0.1;
%	taur_tot = [-0.01, 0.05, dtau];
%	gammar_tot = [0.0, 1.0, dg];
%	taur_tot = [0.008, 0.008, dtau];
%	gammar_tot = [0.7, 0.8, dg];
end
