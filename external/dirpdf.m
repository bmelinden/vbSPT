function p = dirrnd(x,alpha,varargin)
	% p = dirrnd(x,alpha,[N M ...])
	%
	% Computes the proability density function for one or more instances x
	% of a Dirichlet distribution, i.e.,  
	%
	%	p(x|alpha)= beta(k) ~ Gamma(alpha(k), 1)
	%	pi(k) = beta(k) / sum_l beta(l)
	% 
	% If alpha is an N-dimensional array, draws are normalized 
	% along the last dimension.
    %
    % Jan-Willem van de Meent, downloaded 2013-04-12 from 
    % https://raw.github.com/jwvdm/lspace/master/code/matlab/dirrnd.m

	% determine normalization axis
	d = find(size(alpha) > 1, 1, 'last');

	if nargin > 1  
		% replicate alpha if necessary
		dims = [cat(2, varargin{:}) 1];
		o = ones(dims);
		alpha = bsxfun(@times, o, reshape(alpha, [ones(1,ndims(o)) size(alpha)]));
		% adjust normalization dim
		d = d + ndims(o);
	end

	a = gamrnd(alpha, 1);
	p = squeeze(bsxfun(@rdivide, a, sum(a, d)));