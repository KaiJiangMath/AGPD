%% plot phase diagram
function diagram(model, pattern_tot, num)

%model = 'lp';
%pattern_tot = 'lam hex';
%num = '0';

	%clear all; clc;
	close all;
	set(0, 'DefaultFigureVisible', 'on');

	%% the file path
	fstr = sprintf('%s_results/', model);
	doc = strcat(fstr, 'diagram/');
	if exist(doc) == 0
		mkdir(doc);
	end

	%% add the function path
	addpath(genpath('src/'));

	%% set the name of x- and y-axis
	if strcmp(model, 'lp') || strcmp(model, 'multiscales')
		xynotes = {'$\epsilon$', '$\alpha$'};
		type = [2, 3]; % x: 2nd; y: 3rd
	else
		xynotes = {'$\tau$', '$\gamma$'};
		type = [2, 3]; % x: 2nd; y: 3rd
	end

	%% obtain the data of all patterns
	pattern_tot = cellstr(regexp(pattern_tot, ' ', 'split'));   % string to cell
	len = length(pattern_tot);
	data = cell(1, len);
	diagram_data = cell(1, len);
	for i = 1:1:len
		pattern = pattern_tot{i};
		fstr_pattern = sprintf('%s%s/%s_hamilton.txt', fstr, pattern, pattern);
		data{i} = load(fstr_pattern);
	end

	%% deal with the load data (maybe the first, second columns are not equal for all patterns)
	parameters_tot = [];
	for i = 1:1:len
		if ~isempty(data{i})
			parameters_tot = [parameters_tot; data{i}(:,1:1:2)]; % load tau and gamma
		end
	end	
	parameters_tot = mysort(parameters_tot); % sort; first by tau, then by gamma
	%% remove the duplicate data
	coordinates = parameters_tot(1,:);
	for ij = 2:1:size(parameters_tot,1)
		err = parameters_tot(ij-1,:) - parameters_tot(ij,:);
		if ( max(abs(err)) > 1.0e-8 )
			coordinates(end+1,:) = parameters_tot(ij,:);
		end
	end
	[nr, nc] = size(coordinates);

	%% clear up data
	%% in the computed data, the coordinates of different patterns maybe different
	for i = 1:1:len
		temp1 = data{i}; % temporal data for hamilton
		if ~isempty(temp1)
			temp2 = zeros(nr, size(temp1,2));
			for j1 = 1:1:nr
                flag = 1;
				for j2 = 1:1:size(temp1)
					err = coordinates(j1,:) - temp1(j2,1:1:2);
					if ( max(abs(err)) < 1.0e-8 )
                        flag = 0;
						temp2(j1,:) = temp1(j2,:);
						break; % break j2 recurrence
					end
				end
                if ( flag == 1 ) % no data for the set of parameter
                    temp2(j1,:) = [coordinates(j1,:), 10, 0];
                end
			end
		else
			temp2 = [coordinates, 10*ones(nr,1), zeros(nr,1)];
		end
		% the data for identifying diagram
		diagram_data{i} = temp2;
		% deal hamilton according to flag;
		data{i} = dealData(temp2);
%		data{i} = temp2;
	end

	%% save the data of diagram;
	fname = sprintf('%sdiagram%s.mat', doc, num);
	save(fname, 'data');


	%% my defining colors
	mycolor = [
		28, 28, 28; % grey11
		255, 0, 0; % Red1
		0, 139, 139; % DarkCyan
		139, 0, 139; % DarkMagenta
		144, 238, 144; % lightGreen
		139, 0, 0; % DarkRed
		191, 62, 255; % DarkOrchid1
		0, 0, 139; % DarkBlue
		255, 69, 0; % OrangeRed1
		255, 185, 15; % DarkGoldenrod1
		47, 79, 79; % DarkSlateGray
		0, 191, 255; % DeepSkyBlue
		0, 100, 0; % DarkGreen
		255, 255, 0; % Yellow
		178, 34, 34; % Firebrick
		255, 165, 79; % Tan1
		]/255;

	%% adjust the size of marker
	if size(data{1},1) > 3000
		markersize = 10;
	elseif size(data{1},1) > 2000
		markersize = 15;
	elseif size(data{1},1) > 1000
		markersize = 20;
	elseif size(data{1},1) > 500
		markersize = 25;
	else
		markersize = 30;
	end
	markersize = 6;


	%% plot the legend of the phase diagram
	figure('Name', 'legend', 'NumberTitle', 'off', 'Position', [400,20,1200,800]);
	set(gca, 'LineWidth', 2, 'FontName', 'Times New Roman', 'FontSize', 25);
	hold on;
	box off;
	axis off;
	for i = 1:1:size(mycolor,1)
		plot(0, 0, 'LineStyle', 'none', 'Marker', 'o',...
			'MarkerSize', 25, 'MarkerFaceColor',...
			mycolor(i,:), 'MarkerEdgeColor', mycolor(i,:));
	end
	legend(pattern_tot, 'Location', 'north', 'FontName', 'Times New Roman', 'FontSize',...
		30, 'FontWeight', 'bold', 'LineWidth', 2);
	xlabel('');
	ylabel('');
	saveas(gcf, [doc, 'legend.png'], 'png');
	fprintf('finish legend\n');


	%% plot phase diagram
	figure('Name', 'phase diagram', 'NumberTitle', 'off', 'Position', [400,20,1200,800]);
	set(gca, 'LineWidth', 2, 'FontName', 'Times New Roman', 'FontSize', 25);
	hold on;
	box on;
	axis tight;

	for i = 1:1:size(data{1},1)
		energy = zeros(1, len);
		for j = 1:1:len
			energy(j) = data{j}(i,3);
		end
		[val,ind] = min(energy);

		if ( val < -5.0e-8 )
%		if ( val < -1.0e-6 ) % ordered patterns
			diagram_data{ind}(i, 4) = 1; % the stable pattern in this position
			plot(data{1}(i,type(1)-1), data{1}(i,type(2)-1), 'LineStyle', 'none',...
				'Marker', 'o', 'MarkerSize', markersize, 'MarkerFaceColor',...
				mycolor(ind,:), 'MarkerEdgeColor', mycolor(ind,:));
			for j = 1:1:len % others lost
				if ( j ~= ind && diagram_data{j}(i,4) ~= 2 )
					diagram_data{j}(i,4) = 0;
				end
			end
		else % disorder patterns % no order patterns take the dominant role
			for j = 1:1:len
				diagram_data{j}(i,4) = 0;
			end
			plot(data{1}(i,type(1)-1), data{1}(i,type(2)-1), 'LineStyle', 'none',...
				'Marker', 'x', 'MarkerSize', markersize, 'MarkerFaceColor',...
				[0,0,0], 'MarkerEdgeColor', [0,0,0]);
		end
	end

	%% annotation of x and y
	%xlim([0, 2.8])
	%ylim([-0.4,0.3])
	xlabel(xynotes{type(1)-1}, 'Interpreter', 'Latex', 'FontName',...
		'Times New Roman', 'FontSize', 50);
	ylabel(xynotes{type(2)-1}, 'Interpreter', 'Latex', 'FontName',...
		'Times New Roman', 'FontSize', 50, 'Rotation', pi/2);

	figfile = sprintf('%sdiagram%s.png', doc, num);
	saveas(gcf, figfile, 'png');
	fprintf('finish diagram\n');


	%% save the data of the phase diagram
	for j = 1:1:len
		pattern = pattern_tot{j};
		fstr_pattern = sprintf('%s%s/%s_diagram_hamilton.txt', fstr, pattern, pattern);
		eham = fopen(fstr_pattern, 'w');
		for i = 1:1:size(data{j},1)
			fprintf(eham, '% .6f\t% .6f\t ', diagram_data{j}(i,1), diagram_data{j}(i,2));
			fprintf(eham, '% .6e\t%d\n', diagram_data{j}(i,3), diagram_data{j}(i,4));
		end
		fclose(eham);
	end
	fprintf('finish the data of diagram\n');

%
%	%% save the data of the phase boundary;
%	for j = 1:1:len
%		index = find(diagram_data{j}(:,4)==1);
%		coordinate = diagram_data{j}(index, 1:1:2);	% screen the coordinates of j-phase
%		ind = boundary(coordinate(:,1), coordinate(:,2)); % find the index of the boundary of j-phase
%		boundary_coordinate = coordinate(ind, :); % obtain the coordinates of the boundary of j-phase
%		pattern = pattern_tot{j};
%		fstr_pattern = sprintf('%s%s/%s_boundary.txt', fstr, pattern, pattern);
%		eham = fopen(fstr_pattern, 'w');
%		for i = 1:1:size(boundary_coordinate,1)
%			fprintf(eham, '% .6f\t% .6f\n', boundary_coordinate(i,1), boundary_coordinate(i,2));
%		end
%		fclose(eham);
%		fprintf('finish the data of %s phase boundary\n', pattern);
%	end
%
	
	%% remove the function path
	rmpath(genpath('src/'));

	%% the code is finished
	finish = fopen('finish/diagram.txt', 'w');
	fprintf(finish, 'finish');
	fclose(finish);
end
