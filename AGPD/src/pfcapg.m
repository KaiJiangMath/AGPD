%% iteration method adopts the APG algorithm
function [hamilton, uCplx] = pfcapg(model, flow_type, PATTERN,...
	k2, sysPmts, dim, cmpPmts, rcpBox, uc0, fname)

      APGDiffTol=1e-15;
%     format long e
%     if (nargin<7) ||(isempty(APGDiffTol))
%         APGDiffTol = 1e-15;     % initial the APGDiffTol
%     end

  	DimPhy = dim(1);
	DimCpt = dim(2);

%    lid = fopen(fname, 'a');

	%%%%%%% model systems %%%%%%%
	scaleQ = sysPmts(1);   %%%  scale 
	paraQ  = sysPmts(2);   %%%  another scale Q
	paraC  = sysPmts(3);   %%%  interaction term
	paraXi = sysPmts(4);   %%%  2nd term
	paraA  = sysPmts(5);   %%%  3rd term
	paraG  = sysPmts(6);   %%%  4th term
	%%%%%%% model systems %%%%%%%

	dirBox = getDualBox(rcpBox, scaleQ);

	TOL = cmpPmts(1);
	itMax = cmpPmts(2);
	tstep = cmpPmts(3);
	tmin = cmpPmts(4);
	tmax = cmpPmts(5);
	rho = 0.5*(sqrt(5)-1);
	delta = 1e-14;
    
    %Nc = 5;%%% record the # of histories
	t0 = tstep;
	hamilton = inf;
	iterator = 0;
	err = 1.0;

	ncpt = size(uc0);

	q2 = k2;
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
	if strcmp(flow_type, 'AC')
		flow = 1;
	elseif strcmp(flow_type, 'CH')
		flow = q2;
	end

	uCplx = uc0;
    uReal = real(fftn(uCplx));
%%%%%  For Nesterov accelerated method %%%%%%%%
	yCplx = uCplx;
	theta = 1.0;
	q = 0.0;
    RestartFlag = 0;
%%%%%  For Nesterov accelerated method %%%%%%%%
    
    %% Main part %%%%%%%%%%%%%%%%%
	while(abs(err) > TOL && iterator < itMax)
		
		iterator = iterator+1;
		
        %%%%% compute the imformation of extrapolation step %%%%%
		yReal = real(fftn(yCplx)); 
		y2Real = yReal.*yReal; y3Real = y2Real.*yReal; y4Real = y3Real.*yReal;
		gyeReal = paraXi*yReal - 0.5*paraA*y2Real + paraG*y3Real/6;
		gyeCplx = myifftn(gyeReal);

		yeReal = paraXi*y2Real/2 - paraA*y3Real/6 + paraG*y4Real/24;
		yeCplx = myifftn(yeReal);
        intmpyCplx = poten.*(abs(yCplx)).^2;
		intmpyReal = fftn(intmpyCplx);
		yham = real(yeCplx(1)) + (0.5)*real(intmpyReal(1));

		
	%%%%%%  estimate step length  %%%%%%	
	    if iterator == 1 || RestartFlag
            RestartFlag = 0;
			tstep = t0;
		else
			%%%%%  BB step %%%%%
			gyediff = gyeCplx - gyec0;
			yediff  = yCplx - yc0;
            
            %if mod(iterator, 2) == 0
            %   a1 = sum(k2(:).*yediff(:).^2)/sum(k2(:).*yediff(:).*gyediff(:));
            %else
            %   a1 = sum(k2(:).*yediff(:).*gyediff(:))/sum(k2(:).*gyediff(:).^2);
            %end
            
            if mod(iterator, 2) == 0
               a1 = sum(yediff(:).^2)/sum(yediff(:).*gyediff(:));
            else
		   	   a1 = sum(yediff(:).*gyediff(:))/sum(gyediff(:).^2);
            end
            
			tstep = abs(real(a1));
            tstep = min(max(tstep, tmin),tmax);
        end
        
        %%%%%%%%% line search %%%%%%%%%%%%
		while 1 
				zCplx = (yCplx - tstep*flow.*gyeCplx)./(1+tstep*flow.*poten);
                zCplx(1) = 0;
                
                zReal = real(fftn(zCplx)); 
				z2Real = zReal.*zReal; z3Real = z2Real.*zReal; z4Real = z3Real.*zReal;


				zeReal = 0.5*paraXi*z2Real - paraA*z3Real/6 + paraG*z4Real/24;
				zeCplx = myifftn(zeReal);
                
                intmpCplx = poten.*(abs(zCplx)).^2;
                intmpReal = fftn(intmpCplx);
				zham = real(zeCplx(1)) + 0.5*real(intmpReal(1));			

                
				tmp = zReal-yReal;
%                 tmp = zCplx-yCplx;
				tmp1 = yham - zham - delta*norm(tmp(:), 2);
				tmp = uReal-yReal;
%                 tmp = uCplx - zCplx;
				tmp2 = hamilton - zham - delta*norm(tmp(:), 2);
% 				[tmp1, tmp2]
				if (tmp1 <= 0) && (tmp2 <= 0) && (tstep > tmin)
					tstep = rho*tstep;
				else
% 					flag = 0; 
					break; 
                end
        end %%%% end while
   
    %%	
        diffham = zham-hamilton;
        
		if(diffham > 0)
        %%%%%%  restart parameter %%%%%%
			theta = 1;
            yCplx = uCplx;
		% disp('Restart');
            RestartFlag = 1;
        else
        %%%%%%  accept this step %%%%%%
            hamilton = zham;
            ucTmp = uCplx;
            uCplx = zCplx;
            
            gzeReal = paraXi*zReal - 0.5*paraA*z2Real + paraG*z3Real/6;
            gzeCplx = myifftn(gzeReal);  
            gc = poten.*zCplx + gzeCplx;
            gc(1)  = 0;
            err = max(abs(gc(:)));
            
            %%%% update the extrapolation  step %%%%
            theta1 = -0.5*(theta^2-q) + sqrt(0.25*(theta^2-q)^2 + theta^2);
            beta = theta*(1-theta) / (theta^2+theta1);
            theta = theta1;
            
            yc0 = yCplx;
            gyec0 = gyeCplx;
            
            yCplx = (1.0+beta)*uCplx - beta*ucTmp;
        end
        

	%%%%%%%%% Record the information %%%%%%%%	   
		if mod(iterator, 100) == 0
			fprintf('\t --> step %d : error = %.10e\t tstep = %.3e\t hamilton = %.10e \t diffham = %.10e \n', iterator, err, tstep, hamilton, diffham);
		end
	
       if abs(diffham) < APGDiffTol
           break
       end
	end  %%  END while 
%    fclose(lid);

	fprintf('\t --> step %d : error = %.10e\t tstep = %.3e\t hamilton = %.10e \t diffham = %.10e \n', iterator, err, tstep, hamilton, diffham);
end
