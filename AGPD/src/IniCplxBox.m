%% obtain the initial reciprocal value and box;
function [uinit, rcpBox, dim] = IniCplxBox(PATTERN, nc, scaleQ)

	%%%%   dim(1) = DimCpt; dim(2) = DimPhy;

	if strcmp(PATTERN, 'lam')
		dim = [2,2];
		L = 1;
		rcpBox = zeros(dim(1));
		for j=1:1:dim(1) rcpBox(j,j)=scaleQ/L; end
		uinit = obtainIniVals(PATTERN, nc, dim);
	elseif strcmp(PATTERN, 'gyroid')
		dim = [3,3];
		L = sqrt(5);
		rcpBox = zeros(dim(1));
		for j=1:1:dim(1) rcpBox(j,j)=scaleQ/L; end
		uinit = obtainIniVals(PATTERN, nc, dim);
	elseif strcmp(PATTERN, 'hex')
		dim = [2,2];
		rcpBox = [1, cos(pi/3);
				0, sin(pi/3);];
		uinit = obtainIniVals(PATTERN, nc, dim);
	elseif strcmp(PATTERN, 'hcp')
		dim = [3,3];
		L = sqrt(2);
		rcpBox = zeros(dim(1));
		for j=1:1:dim(1) rcpBox(j,j)=scaleQ/L; end
		uinit = obtainIniVals(PATTERN, nc, dim);
	elseif strcmp(PATTERN, 'bcc')
		dim = [3,3];
		L = sqrt(2);
		rcpBox = zeros(dim(1));
		for j=1:1:dim(1) rcpBox(j,j)=scaleQ/L; end
		uinit = obtainIniVals(PATTERN, nc, dim);
	elseif strcmp(PATTERN, 'fcc')
		dim = [3,3];
		L = sqrt(3);
		rcpBox = zeros(dim(1));
		for j=1:1:dim(1) rcpBox(j,j)=scaleQ/L; end
		uinit = obtainIniVals(PATTERN, nc, dim);
	elseif strcmp(PATTERN, 'A15')
		dim = [3,3];
		L = sqrt(5);
		rcpBox = zeros(dim(1));
		for j=1:1:dim(1) rcpBox(j,j)=scaleQ/L; end
		uinit = obtainIniVals(PATTERN, nc, dim);
	elseif strcmp(PATTERN, 'fddd')
		dim = [3,3];
		dirBox = [  pi, 0, 0; 
					0, 2*pi, 0; 
					0, 0, 2*sqrt(3)*pi
				];
		rcpBox = getDualBox(dirBox, scaleQ);
		uinit = obtainIniVals(PATTERN, nc, dim);
	elseif strcmp(PATTERN, 'fdddCube')
		dim = [3,3];
		L = 1/7.92;
		rcpBox = zeros(dim(1));
		for j=1:1:dim(1) rcpBox(j,j)=scaleQ/L; end
		uinit = obtainIniVals(PATTERN, nc, dim);
	%%   load from file
	%    nc = [24,24,24];
	%    tmp = load('./initvals/FdddCube/fdddcubeA');
	%    ur0 = reshape(tmp(4:end,:), nc);
	%    uinit = myifftn(ur0);
	%
	%    sflag = detsymmetry(uinit, uinit);
	%
	%	 dirBox = diag([tmp(end, 1), tmp(end,2), tmp(end,3)]);
	%    rcpBox = getDualBox(dirBox, scaleQ);
	%    clear tmp;
	elseif strcmp(PATTERN, 'sigma')
		dim = [3,3];
		if nc(end) == 64
			tmp = load('./initvals/SigmaPhase/Sigma128x128x64.txt');
		elseif nc(end) == 128
			tmp = load('./initvals/SigmaPhase/Sigma256x256x128.txt');
		end
		ur0 = reshape(tmp(:,4), nc);
		uinit = myifftn(ur0);
		%%
		dirBox = diag([tmp(end, 1), tmp(end,2), tmp(end,3)]);
		rcpBox = getDualBox(dirBox, scaleQ);
		clear tmp;
	elseif strcmp(PATTERN, '10fold')
		dim = [2,4];
		L = 1;
		rcpBox = zeros(dim(1));
		for j=1:1:dim(1) rcpBox(j,j)=scaleQ/L; end
		uinit = obtainIniVals(PATTERN, nc, dim);
	elseif strcmp(PATTERN,'12fold') || strcmp(PATTERN, 'LQ6') || strcmp(PATTERN, 'LQS6') ||...
			strcmp(PATTERN, 'C3') || strcmp(PATTERN, '12i6o') || strcmp(PATTERN, '8i10o') ||...
		   	strcmp(PATTERN, 'sq') || strcmp(PATTERN, 'squ') || strcmp(PATTERN, 'sqv') ||...
			strcmp(PATTERN, 'sqw') || strcmp(PATTERN, 'Ls')
		dim = [2,4];
		L = 1;
		rcpBox = zeros(dim(2));
		for j=1:1:dim(2) rcpBox(j,j)=scaleQ/L; end
		uinit = obtainIniVals(PATTERN, nc, dim);
	end
end
