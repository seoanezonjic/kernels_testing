function msg = kernes(kernel_type,matrix_file,outf)
	% Load matrix
	data = matfile(matrix_file)
	A = data.M
	D = sum(A,2);

	% Perform target kernel
	if strcmp('ct',kernel_type)
		L = diag(D)-A;
		M = pinv(L);
		save(outf,'M');
		msg = 'Commute-Time kernel ended without errors'
	elseif strcmp('el',kernel_type)
		beta=0.02;
		L = diag(D)-A;
		M = expm(-beta*L);
		save(outf,'M');
		msg = 'Laplacian Diffusion kernel ended without errors'
	elseif strcmp('ka',kernel_type)
		lambda = min(eig(A));
		M = A + eye(length(D))*abs(lambda);	
		save(outf,'M');
		msg = 'Kernelized Adjacency matrix ended without errors'
	elseif strcmp('me',kernel_type)
		beta=0.04;
		T = (eye(length(D))*length(D) - diag(D) + A) * (beta/length(D));
		M = expm(T);
		save(outf,'M');
		msg = 'Markov exponential diffusion kernel ended without errors'
	elseif strcmp('rf',kernel_type)
		L = diag(D)-A;
		M = inv(eye(length(D)) + L)
		save(outf,'M');
		msg = 'Random forest kernel ended without errors'
	elseif regexp(kernel_type,'vn*')
		p = str2num(regexprep(kernel_type,'vn',''))
		alpha = p * (max(eig(A))^(-1));
		M = inv(eye(length(D)) - A*alpha);
		save(outf,'M');
		msg = 'Von-Neumann kernel ended without errors'
	else
		msg = 'Kernel type not supported'
	end
end