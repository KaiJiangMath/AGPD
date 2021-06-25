%% identify pattern by perception Hash algorithm
%% return a Hash vector according to the picture fname
function hash_feature = aHash(fname, N)

	I = imread(fname);
	%% scaling picture
	J = imresize(I, [N,N]);
	%% N2 gray, coarse
	img = double(rgb2gray(J)); 
	img_N2 = floor(img/255*(N^2));
	%% the mean of the gray values
	gray_mean = sum(img_N2(:))/(N^2);
	%% compare the gray values of this picture with the mean value
	%% 1: greater or equal; 0: smaller
	feature_img = zeros(N, N);
	for i = 1:1:N
		for j = 1:1:N
			if img_N2(i,j) >= gray_mean
				feature_img(i,j) = 1;            
			end
		end
	end
	%% calculate the Hash value
	%% compare the Hash value to discuss the similarity of pictures
	hash_feature = reshape(feature_img, 1, N^2);
end
