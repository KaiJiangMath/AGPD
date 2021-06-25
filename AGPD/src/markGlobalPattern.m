%% obtain coordinates and values in Fourier space
function kindex = markGlobalPattern(PATTERN)

	if strcmp(PATTERN, 'lam')
		kindex = [ 
			1		0		0.058;
			-1		0		0.058;
			];
	elseif strcmp(PATTERN, '2fold')
		kindex = [
			1,		0.01; 
			-1,		0.01; 
			];
	elseif strcmp(PATTERN, 'hex')
		kindex = [
				1     0     0.3;
			   -1     0     0.3;
				0     1     0.3;
				0    -1     0.3;
				1    -1     0.3;
			   -1     1     0.3;
			];
	elseif strcmp(PATTERN, 'hcp')
		kindex = [
				0     1    -1	0.3; 
				0    -1     1	0.3; 
				1     0    -1	0.3; 
				1    -1     0	0.3; 
			   -1     0     1	0.3; 
			   -1     1     0	0.3; 
			];
	elseif strcmp(PATTERN, 'LQ6') || strcmp(PATTERN, 'LQS6') || strcmp(PATTERN, 'C3') ||...
	   	strcmp(PATTERN, '12i6o') || strcmp(PATTERN, '8i10o') || strcmp(PATTERN, 'sq') ||...
		strcmp(PATTERN, 'squ') || strcmp(PATTERN, 'sqv') || strcmp(PATTERN, 'sqw') ||...
		strcmp(PATTERN, 'Ls')
		kindex = [
			 1,  0,  0,  0;
			 0,  1,  0,  0;
			 0,  0,  1,  0;
			 0,  0,  0,  1;
			-1,  0,  1,  0;
			 0, -1,  0,  1;
			 1,  1,  0,  0;
			 0,  1,  1,  0;
			 0,  0,  1,  1;
			-1,  0,  1,  1;
			-1, -1,  1,  1;
			-1, -1,  0,  1;
			];
		switch PATTERN
		case 'LQ6'
			index = [2:2:6, 1, 7, 12];				%% 4 LQ6;
		case 'LQS6'
			index = [8:2:12, 11, 2, 3];				%% 5 LQS6 7:2:11 10 1 2
		case 'C3'
			index = [2, 3, 8];						%% 6 C3;
		case '12i6o'
			index = [1:6, 7:2:11];					%% 9 12i6o;
		case '8i10o'
			index = [2, 3, 5, 6, 7, 9, 10, 11, 12];	%% 10 8i10o;
		case 'sq'
			index = [2, 3, 5, 6, 8, 11];			%% 11 sq;
		case 'squ'
			index = [1, 4];							%% 12 squ;
		case 'sqv'
		    index = [1, 4, 7, 12];					%% 13 sqv;
		case 'sqw'
			index = [1, 4, 6, 7];					%% 14 sqw;
		case 'Ls'
			index = [7];							%% 15 Ls;
		end
		kindex = kindex(index, :);
		kindex = [kindex; (-kindex)];
		kindex = [kindex, 0.058 * ones(size(kindex,1),1)];
	elseif strcmp(PATTERN, '12fold')
		kindex = [
			0     1     0   -1	0.058;
			0    -1     0    1	0.058;
			1     0     0    0	0.058;
		   -1     0     0    0	0.058;
			0     1     0    0	0.058;
			0    -1     0    0	0.058;
			0     0     1    0	0.058;
			0     0    -1    0	0.058;
			0     0     0    1	0.058;
			0     0     0   -1	0.058;
		   -1     0     1    0	0.058;
			1     0    -1    0	0.058;
			1     1     0   -1	0.058;
		   -1    -1     0    1	0.058;
			1     1     0    0	0.058;
		   -1    -1     0    0	0.058;
			0     1     1    0	0.058;
			0    -1    -1    0	0.058;
			0     0     1    1	0.058;
			0     0    -1   -1	0.058;
		   -1     0     1    1	0.058;
			1     0    -1   -1	0.058;
		   -1    -1     1    1	0.058;
			1     1    -1   -1	0.058;
			];
	elseif strcmp(PATTERN, '10fold')
		kindex = [      
			1	  0     0    0	0.058;
			0	  1     0    0	0.058;
			0	  0     1    0	0.058;
			0	  0     0    1	0.058;
			-1	  1    -1    1	0.058;
			-1	  0     0    0	0.058;
			0	 -1     0    0	0.058;
			0	  0    -1    0	0.058;
			0	  0     0   -1	0.058;
			1	 -1     1   -1	0.058;
			1	  0     1   -1	0.058;
			1	  0     1    0	0.058;
			0	  1     0    1	0.058;
			-1	  1     0    1	0.058;
			-1	  0     0    1	0.058;
			-1	  0    -1    1	0.058;
			-1	  0    -1    0	0.058;
			0	 -1     0   -1	0.058;
			1	 -1     0   -1	0.058;
			1	  0     0   -1	0.058;
			];
	elseif strcmp(PATTERN, 'CAM12fold')
		kindex = [
			30     0	0.058;
			26    15	0.058;
			15    26	0.058;
			0     30	0.058;
		   -15    26	0.058;
		   -26    15	0.058;
		   -30     0	0.058;
		   -26   -15	0.058;
		   -15   -26	0.058;
			0    -30	0.058;
			15   -26	0.058;
			26   -15	0.058;
			56    15	0.058;
			41    41	0.058;
			15    56	0.058;
		   -15    56	0.058;
		   -41    41	0.058;
		   -56    15	0.058;
		   -56   -15	0.058;
		   -41   -41	0.058;
		   -15   -56	0.058;
			15   -56	0.058;
			41   -41	0.058;
			56   -15	0.058;
			];
	elseif strcmp(PATTERN, 'cam6fold')
		kindex = [
			30     0	0.1;
			15    26	0.1;
		   -15    26	0.1;
		   -30     0	0.1;
		   -15   -26	0.1;
			15   -26	0.1;
			];
	elseif strcmp(PATTERN, 'bcc')
		kindex = [
			1	 1	 0	0.3;
		   -1	 1	 0	0.3;
		   -1	-1	 0	0.3;
			1	-1	 0	0.3;
			0	 1	 1	0.3;
			0	-1	 1	0.3;
			0	-1	-1	0.3;
			0	 1	-1	0.3;
			1	 0	 1	0.3;
		   -1	 0	 1	0.3;
		   -1	 0	-1	0.3;
			1	 0	-1	0.3;
			];
	elseif strcmp(PATTERN, 'fcc')
		kindex = [
			1	 1	 1	0.3;
			1	-1	 1	0.3;
		   -1	 1	 1	0.3;
		   -1	-1	 1	0.3;
			];
	elseif strcmp(PATTERN, 'A15')
		kindex = [
		   -2, -1, 0 	0.0338277052131978;
			2,  1, 0 	0.0338277052131978;
			2, -1, 0 	0.0338277052131977;
		   -2,  1, 0 	0.0338277052131977;
		   -1,  2, 0   -0.0338277052131977;
			1, -2, 0   -0.0338277052131977;
			1,  2, 0   -0.0338277052131977;
		   -1, -2, 0   -0.0338277052131977;
			0, -2, 1 	0.0338228603092057;
		   -2,  0, 1   -0.0338228603092056;
			0,  2, 1 	0.0338228603092055;
			2,  0, 1   -0.0338228603092055;
		   -1,  0, 2 	0.0338139694707322;
			0, -1, 2   -0.0338139694707320;
			0,  1, 2   -0.0338139694707319;
			1,  0, 2 	0.0338139694707319;
			0,  0, 2 	0.0331853928814813;
			0,  2, 0 	0.0331847214985317;
			2,  0, 0 	0.0331847214985317;
			0, -2, 0 	0.0331847214985317;
		   -2,  0, 0 	0.0331847214985317;
			1, -1, 2 	0.0297069055330425;
		   -1, -1, 2 	0.0297069055330425;
			1,  1, 2 	0.0297069055330424;
		   -1,  1, 2 	0.0297069055330424;
		   -1,  2, 1 	0.0297065456817734;
		   -2, -1, 1 	0.0297065456817733;
		   -1, -2, 1 	0.0297065456817733;
			2, -1, 1 	0.0297065456817733;
		   -2,  1, 1 	0.0297065456817732;
			2,  1, 1 	0.0297065456817732;
			1,  2, 1 	0.0297065456817732;
			1, -2, 1 	0.0297065456817731;
			];
	elseif strcmp(PATTERN, 'gyroid')
		kindex = [
			 1	-2	1	 7.433127e-02;
			-2	 1	1	 7.433126e-02;
			-1	-1	2	 7.433125e-02;
			 1	 2	1	-7.433124e-02;
			 1	 1	2	-7.433123e-02;
			 2	 1	1	-7.433122e-02;
			-2	-1	1	-7.433121e-02;
			 1	-1	2	-7.433120e-02;
			-1	 1	2	 7.433120e-02;
			-1	 2	1	-7.433119e-02;
			-1	-2	1	 7.433119e-02;
			 2	-1	1	 7.433119e-02;
			];
	elseif strcmp(PATTERN, 'fddd')
		kindex = [
				1     1     1  -0.3; 
				1     1    -1	0.3; 
				1    -1     1	0.3; 
				1    -1    -1	0.3;
				0     2     2  -0.3;
				0     2    -2   0.3;
				0     0     4	0.3;
			];
	elseif strcmp(PATTERN, 'fdddCube')
		kindex = [
				2	-1	 1   0.06;  
				1	-1	 2  -0.06;  
			   -1	 1	-2  -0.06;  
			   -2	 1	-1   0.06;  
				2	 1	 1  -0.06;  
				1	 1	 2  -0.06;  
			   -1	-1	-2  -0.06;  
			   -2	-1	-1  -0.06;  
			   -1	 2	 1  -0.03;  
			   -1	-2	 1   0.03;  
				1	 2	-1   0.03;  
				1	-2	-1  -0.03;  
			   -2	 0	 2   0.02;  
				2	 0	-2   0.02;  
				0	 1	 3  -0.005;  
				0	-1	 3   0.005;  
				0	 1	-3   0.005;  
				0	-1	-3  -0.005;  
				3	 1	 0  -0.005;  
			   -3	 1	 0  -0.005;  
				3	-1	 0  -0.005;  
			   -3	-1	 0  -0.005
		];
	elseif strcmp(PATTERN,'icosahedron')
		tmp = [
				1	0	0	0	0	0; 
				0	1	0	0	0	0; 
				0	0	1	0	0	0; 
				0	0	0	1	0	0; 
				0	0	0	0	1	0; 
				0	0	0	0	0	1; 
				0	0	1	0	1	0;
				0	1	0	0  -1	0;
				0	1  -1	1  -1	1;
				1  -1   1  -1   1  -1;
				1	0	0  -1	0	0; 
			   -1	1	0   1	0	1; 
			   -1	1	0   1	0	0; 
				0	0   1  -1	0  -1; 
				0	0   1  -1	0   0;
				]; 
		kindex = [tmp; -1.0*tmp];
		tmpq = [
			   -1	1	1	0	0	0;
				0	0	0	1	0	1;
				0	1  -1	1  -1	0;
				1	0	0  -1	0  -1;
				1  -1	1  -1	1	0;
				0	0	1  -1	1	0;
				0	0   1  -1	1  -1;
			   -1   1	0	1  -1	0;
			   -1	1	0	1  -1	1;
				1	0  -1	0	0	0;
				1  -1	0	0	1	0;
				1  -1	0	0	0	0;
				1	0  -1	0  -1	0;
			   -1	1  -1	2  -1	1;
				0	1  -1	0  -1	0;
				];
		kindex = [kindex; tmpq; -1.0*tmpq];
		val = 0.03 * ones(size(kindex,1), 1);
		kindex = [kindex, val];
	else
		fidtxt = sprintf('initvals/nfoldQC/%s.txt', PATTERN);
		kindex = textread(fidtxt);
		kindex = [kindex, 0.058*ones(size(kindex,1),1)];
	end
end
