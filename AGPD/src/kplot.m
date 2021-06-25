
function kplot(kk, dim)

DimPhy = dim(1);

figure(1); hold on;
axis off
axis equal;
set(gcf, 'color', 'white')
colormap;
set(gcf, 'unit', 'normalized', 'position', [0.05, 0.1, 0.6, 0.6])
title('Spectral Distribution of Density')
xlabel('')
ylabel('')

%kspace = sortrows(kk, DimPhy+1, 'descend');
[val, ind] = sort(abs(kk(:,end)), 'descend');
kspace = kk(ind, :);

if DimPhy == 2

	n = length(kspace(:,1));
	x = kspace(:,1);
	y = kspace(:,2);
	index = abs(kspace(:,3));
	plot(0, 0, 'MarkerFaceColor',[0 0 1], ...
	'MarkerEdgeColor',[0 0 1],'Marker','o', ...
	'MarkerSize', 12)
	for i = 1:1:12*8
		if index(i) > 1.0e-2
			plot(x(i),y(i),'MarkerFaceColor',[1 0 0], ...
			'MarkerEdgeColor',[1 0 0],'Marker','o', ...
			'MarkerSize', 10)
		elseif (index(i) < 1.0e-2 && index(i) > 1.0e-3)
			plot(x(i),y(i),'MarkerFaceColor',[0 0 1], ...
			'MarkerEdgeColor',[0 0 1],'Marker','o', ...
			'MarkerSize', 8)
		elseif (index(i) < 1.0e-3 && index(i) > 1.0e-4)
			plot(x(i),y(i),'MarkerFaceColor',[0 1 0], ...
			'MarkerEdgeColor',[0 1 0],'Marker','o', ...
			'MarkerSize', 6.5)
		elseif (index(i) < 1.0e-4 && index(i) > 5.0e-5)
			plot(x(i),y(i),'MarkerFaceColor',[112/255 128/255 105/255], ...
			'MarkerEdgeColor',[112/255 128/255 105/255],'Marker','o', ...
			'MarkerSize', 5.5)
		elseif (index(i) < 1.0e-5 && index(i) > 1.0e-8)
			plot(x(i),y(i),'MarkerFaceColor',[0 0 0], ...
			'MarkerEdgeColor',[0 0 0],'Marker','o', ...
			'MarkerSize', 3)
		end
	end
	drawnow;

end

end

