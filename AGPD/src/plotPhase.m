function plotPhase(PATTERN, uc0, reg, dim, rcpBox, fname)

	DimPhy = dim(1);
	DimCpt = dim(2);
	n = size(uc0);

	ncpt = 3*size(uc0);
	uc_rf = zeros(ncpt);
	uc_rf = mymap(uc0, uc_rf, 'cs2rf');

	ur0 = real(fftn(uc_rf));
	fig1 =figure(1);

	if DimCpt == 1 && DimPhy == 1
		plot(real(ur0));
	elseif DimCpt == 2 && DimPhy == 2
		subplot(1,2,1)

		dirBox = getDualBox(rcpBox, 1);   %%% 1 is scaleQ.  
		angle = rad2deg(acos((dirBox(1,1)*dirBox(1,2)+...
				dirBox(2,1)*dirBox(2,2))/(norm(dirBox(:,1))*norm(dirBox(:,2)))));
				theta=deg2rad(angle);
		theta=deg2rad(angle);

		N_0 = length(ur0(:,1));
		N_1=N_0^2;

		D_1=norm(dirBox(:,1));
		D_2=norm(dirBox(:,2));
		ur0 = ur0(:);

		x=linspace(0,D_1-D_1/N_0,N_0);
		y=linspace(0,D_2-D_2/N_0,N_0);
		%% 
		x_0=zeros(1,N_1);
		y_0=zeros(1,N_1);
		%%  affine transformation
		for i=1:N_1
			x_0(i)=x(mod((i-1),N_0)+1)+y(floor((i-1)/N_0)+1)*cos(theta);
			y_0(i)=y(floor((i-1)/N_0)+1)*sin(theta);
		end
		%% plot
		scatter(x_0,y_0, ncpt(1), ur0,'filled','square');

%        imagesc(real(ur0))
		colormap jet;
		colorbar;
		axis square;
		box on;
		axis on;
		set(gcf, 'color', 'white')

		subplot(1,2,2)
		set(gcf, 'color', 'white')
		kspace = projPlane(PATTERN, uc0, rcpBox, n, dim);
		kplot(kspace, dim);
	elseif DimCpt == 3 && DimPhy == 3
		set(gcf, 'color', 'white');
%        ur0 = -ur0;
		minu = min(ur0(:));
		maxu = max(ur0(:));
		isoA = minu+0.7*(maxu-minu);

		x=linspace(0, reg(1,1), ncpt(1));
		y=linspace(0, reg(2,2), ncpt(2));
		z=linspace(0, reg(3,3), ncpt(3));
		[X Y Z] = meshgrid(x, y, z);

		size([X,Y,Z]);

		alpA = 1;
		patch(isosurface(X,Y,Z, ur0, isoA), ...
		'facecolor',[46,169,223]/255, 'FaceAlpha', alpA, 'edgecolor','none');
		patch(   isocaps(X,Y,Z, ur0, isoA, 'enclose'), ...
		'facecolor','none','FaceAlpha', alpA, 'edgecolor','none');
		daspect([1, 1, 1]);
		camup([1, 0, 0]);
		campos([25, -55, 5]);
		axis square;
		axis equal;
		camlight;
		lighting phong;
		box off;
		axis off;

		if strcmp(PATTERN, 'sigma') view(0, 90);
		else view(112, 35);
		end
		box on;
		axis on;
	elseif DimCpt == 4 && DimPhy == 2
		set(gcf, 'color', 'white')
		kspace = projPlane(PATTERN, uc0, rcpBox, n, dim);
		kplot(kspace, dim);
	end

%	set(gcf,'PaperPositionMode','auto'); % make sure that the saved pic are the same as the display fig
%	print(gcf, '-dpng', '-r280', [fname, '.png']);
    saveas(gcf, [fname, '.png'], 'png')
end
