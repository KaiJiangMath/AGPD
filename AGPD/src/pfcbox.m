function [newBox, q2] = pfcbox(model, PATTERN, sysPmts, dim, uCplx, rcpBox, tol, type)
%%%%%  use BB method to optimize periodicity %%%%%

	%% compute the computational box in direct space
	scaleQ = sysPmts(1);   %%%  scale 
	dirBox = getDualBox(rcpBox, scaleQ)

	box = rcpBox;
	box0 = rcpBox;
	DimCpt = dim(1); 

	lbox = zeros(size(rcpBox));
	rbox = zeros(size(rcpBox));

	%%%% the step length for computing difference quotient
	dh = 1.0e-4;
	dbox = zeros(size(rcpBox));
	gradB = zeros(size(rcpBox));

	%% set the stepping box
	if strcmp(type, 'cube') || ...
	   strcmp(type, 'rectangular')
		for i=1:1:length(dbox(:,1))
			dbox(i,i) = dh;
		end
	elseif strcmp(type, 'hex')
		for i=1:1:length(dbox(:,1))
			dbox(i,i) = dh;
		end
		dbox(2,1) = dh;
	end

	iter = 0;
	eps = 1.0;
	t0 = 0.1;

	%% recurrence for updating the computational box
	while(eps > tol && iter < 50)

		iter = iter + 1;

		%% compute the gradient box
		if strcmp(type, 'cube') || strcmp(type, 'rectangular')
			for i = 1:1:DimCpt
				rbox = box;
				rbox(i,i) = rbox(i,i) + dh;
				[fr,q2] = pfcHam(model, PATTERN, sysPmts, dim, rbox, uCplx);
				lbox = box;
				lbox(i,i) = lbox(i,i) - dh;
				[fl,q2] = pfcHam(model, PATTERN, sysPmts, dim, lbox, uCplx);
				gradB(i,i) = (fr-fl) / (2.0*dh);
			end
		elseif strcmp(type, 'hex')
			for i = 1:1:DimCpt
				for j = 1:1:i
					rbox = box;
					rbox(i,j) = rbox(i,j) + dh;
					[fr,q2] = pfcHam(model, PATTERN, sysPmts, dim, rbox, uCplx);
					lbox = box;
					lbox(i,j) = lbox(i,j) - dh;
					[fl,q2] = pfcHam(model, PATTERN, sysPmts, dim, lbox, uCplx);
					gradB(i,j) = (fr-fl) / (2.0*dh);
				end
			end
		end

		if iter == 1
			tstep = t0;
		else
			gdiff = gradB - gradBold;
			udiff = box1 - box0;
			a1 = sum(udiff(:).^2)/sum(udiff(:).*gdiff(:));
			a2 = sum(udiff(:).*gdiff(:))/sum(gdiff(:).^2);
			tstep = a1;
		end
		
		box = box - tstep*gradB;
		
		if iter == 1
			box1 = box; 
		else 
			box0 = box1;
			box1 = box;
		end
		gradBold = gradB;

		eps = max(abs(gradB(:)));

		fprintf('\t\t ===>  iter : %d\t eps = %.15e\t', iter, eps);
		fprintf('hamilton (+dbox): %.12e\t hamilton (-dbox): %.12e\n', fr, fl);
		%% compute the computational box in direct space
		dirBox = getDualBox(rcpBox, scaleQ)
	end

	newBox = box;
end


%% define the function to calculate the hamilton energy
function [hamilton, q2] = pfcHam(model, PATTERN, sysPmts, dim, rbox, uCplx);

	%%%%%%% model systems %%%%%%%
	scaleQ = sysPmts(1);   %%%  scale 
	paraQ  = sysPmts(2);   %%%  another scale Q
	paraC  = sysPmts(3);   %%%  interaction term
	paraXi = sysPmts(4);   %%%  2nd term
	paraA  = sysPmts(5);   %%%  3rd term
	paraG  = sysPmts(6);   %%%  4th term
	%%%%%%% model systems %%%%%%%

	ncpt = size(uCplx);

	projmat = getprojmat(PATTERN, dim(1), dim(2));
	q2 = obtGsquare(ncpt, rbox, dim, projmat);
	if strcmp(model, 'lb')
		q2(1) = 1;		poten = (q2-1).^2;			poten(1) = 0;
	elseif strcmp(model, 'lp')
		q2(1) = 1;		poten = (q2-1).^2.*(q2-paraQ.^2).^2;		poten(1) = 0;
	elseif strcmp(model, 'ok')
		q2(1) = 1;		poten = (q2-1).^2./q2;		poten(1) = 0;
	elseif strcmp(model, 'leibler')
		poten = leibler_potential(q2);
	end
	poten = paraC * poten;

	uReal = real(fftn(uCplx));
	u2Real = uReal.*uReal; u3Real = u2Real.*uReal; u4Real = u3Real.*uReal;

	enReal = 0.5*paraXi*u2Real - paraA*u3Real/6 + paraG*u4Real/24;
	enCplx = myifftn(enReal);
	intmpCplx = poten.*(abs(uCplx)).^2;
	intmpReal = fftn(intmpCplx);

	hamilton = real(enCplx(1)) + 0.5*real(intmpReal(1));
end
