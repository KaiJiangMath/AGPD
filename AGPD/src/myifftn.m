%%%% use fftn+conj to realize ifftn 
%%%% it costs less cpu time than directly using ifftn

function rslt = myifftn(src)
	s = size(src);
	dof = prod(s);
	rslt = conj(fftn(conj(src))/dof);
end
