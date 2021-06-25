function hamilton = getPartialHamilton(PATTERN, sysPmts, dim, box, uCplx)

	%%%%%%% model systems %%%%%%%
	scaleQ = sysPmts(1);   %%%  scale 
	paraC  = sysPmts(2);  %%%  interaction term
	paraXi = sysPmts(3);   %%%  2nd term
	paraA  = sysPmts(4);   %%%  3rd term
	paraG  = sysPmts(5);   %%%  4th term
	%%%%%%% model systems %%%%%%%

	ncpt = size(uCplx);
	potenCplx = obtainPotential(PATTERN, ncpt, box, dim);
    poten2Cplx = potenCplx.^2;

	uReal = real(fftn(uCplx));
	u2Real = uReal.*uReal; u3Real = u2Real.*uReal; u4Real = u3Real.*uReal;
	entropyReal = paraXi*u2Real/2 - paraA*u3Real/6 + u4Real/24;
	entropyCplx = myifftn(entropyReal);
	entropy = real(entropyCplx(1));
	
	intmpCplx = potenCplx.*uCplx;
	intmpReal = fftn(intmpCplx);
	enIntReal = intmpReal.^2;
	enIntCplx = myifftn(enIntReal);
	interaction = 0.5*paraC*real(enIntCplx(1));

	hamilton = entropy+interaction;

end
