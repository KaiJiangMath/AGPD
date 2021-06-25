
function dbox = getDualBox(box, scale)

	dim = length(box(:,1));
	dbox = zeros(dim, dim);
	
	dbox = (2*pi/scale)*inv(box);
%
%    for i = 1:1:dim
%        dbox(i,i) = (2*pi/scale) / box(i,i) ;
%    end

end

