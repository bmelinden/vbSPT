function x = dirrnd(alpha, varargin)
	% x = dirrnd(alpha, [N M ...])
	%
	% Generates one or more draws from a Dirichlet distribution,
	% i.e. a normalized set of draws from the Gamma distribution: 
	%
	%	beta(k) ~ Gamma(alpha(k), 1)
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
	x = squeeze(bsxfun(@rdivide, a, sum(a, d)));