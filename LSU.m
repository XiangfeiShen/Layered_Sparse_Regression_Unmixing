function [X,sre] = LSU(Y, A, opts)

%--------------------------------------------------------------
% Usage:
%
% X = LSU(Y, A, opts)
%
% ------- Input variables -------------------------------------------
%
%  Y - hyperspectral data matrix with dimensions L(bands) x K(pixels)
%
%  A - spectral library matrix with dimentsions L(bands) x m(spectra)
%
%  parameter.
%            * opts.alpha - regularization parameter of sparsity 
%            * opts.beta - regularization parameter of the ASC constraint  
%            * opts.mu - initial augmented Lagrangian penalty parameter, 
%              default: 0.1
%            * opts.maxiter - maximum iteration number, default: 500
%            * opts.imgsize - the image size (row and column) of input
%            bands
%            * opts.xt - the reference gt
%
%  NOTE: PARAMETERS WITH SYMBOL * ARE NECESSARY
%
% ------- Output variables -------------------------------------------
%
% X - the estimated abundance matrix with dimensions m(spectra) x K(pixels)
% sre - SRE values at each iteration
% ---------------------------------------------------------------------
%
% Please see [1] for more details.
%
% Please contact Xiangfei Shen (xfshen95@163.com) to report bugs or
% provide suggestions and discussions for the codes.
%
% ---------------------------------------------------------------------
% version: 1.0 (8-Sep-2023)
% ---------------------------------------------------------------------
%
% Copyright (Sep, 2023):       Xiangfei Shen (xfshen95@163.com/xfshen95@outlook.com)
%                              Xichuan Zhou* (zxc@cqu.edu.cn)
%
% LSU is distributed under the terms of
% the GNU General Public License 2.0.
%
% Permission to use, copy, modify, and distribute this software for
% any purpose without fee is hereby granted, provided that this entire
% notice is included in all copies of any software which is or includes
% a copy or modification of this software and in all copies of the
% supporting documentation for such software.
%
% This software is being provided "as is", without any express or
% implied warranty.  In particular, the authors do not make any
% representation or warranty of any kind concerning the merchantability
% of this software or its fitness for any particular purpose."
%
% Please cite the following paper if this demo helps your research work.
%
% [1] Xiangfei Shen, Lihui Chen, Haijun Liu Member, IEEE, Xi Su, Wenjia Wei, Xia Zhu,
% and Xichuan Zhou* Senior Member, IEEE, "Efficient Hyperspectral Sparse Regression Unmixing with Multilayers",
% IEEE Transactions on Geoscience and Remote Sensing, Early Access. DOI:10.1109/TGRS.2023.3311642
%
%--------------------------------------------------------------

Y=gpuArray(Y);
A=gpuArray(A);
[L, K] = size(Y);
m = size(A, 2);
if size(A, 1) ~= L
    error(['The sizes of hyperspectral data matrix Y and spectral ',...
        'library matrix A are inconsistent!']);
end

if isfield(opts, 'maxiter')
    MaxIter = opts.maxiter;
else
    MaxIter = 500;
end

if isfield(opts, 'epsilon')
    epsilon = opts.epsilon;
else
    epsilon = 1e-5;
end

if isfield(opts, 'xt')
    XT = opts.xt;
end

if isfield(opts, 'mu')
    mu = opts.mu;
else
    mu = 0.1;
end

if isfield(opts, 'alpha')
    alpha = opts.alpha;
else
    error('The parameter alpha is missing!');
end


if isfield(opts, 'imgsize')
    imgsize = opts.imgsize;
else
    error('The parameter imgsize is missing!');
end

if isfield(opts, 'beta')
    beta = opts.beta;
else
    error('The parameter beta is missing!');
end

if isfield(opts, 'verbose')
    verbose = opts.verbose;
else
    verbose = 1;
end

if verbose==1
    figure
end

%%
%---------------------------------------------
%  Initializations
%---------------------------------------------

I=gpuArray(eye(m));
IF = (A'*A + 3*I)^-1;
Omm=gpuArray(ones(m,m));
Omk=gpuArray(ones(m,K));

U = IF*A'*Y;

V1 = A*U;
V2 = U;
V3 = U;
V4 = U;

D1 = V1*0;
D2 = V2*0;
D3 = V3*0;
D4 = V4*0;

%current iteration number
i = 1;

%primal residual 
res_p = inf;

%dual residual
res_d = inf;

%error tolerance
tol = sqrt((3*m + L)/2*K/2)*epsilon;


%%
%---------------------------------------------
%  ADMM iterations
%---------------------------------------------
while (i <= MaxIter) && ((abs(res_p) > tol) || (abs(res_d) > tol))
    if mod(i, 10) == 1
        V10 = V1;
        V20 = V2;
        V30 = V3;
        V40 = V4;
    end
    %update U and V
    U = IF*(A'*(V1 + D1) + (V2 + D2) + (V3 + D3) + (V4 + D4));
    V1 = 1/(1+mu)*(Y + mu*(A*U - D1));
    V2=solveV2(U-D2,alpha/mu,imgsize);
    V3 =  max(U - D3, 0);%%
    IF2= (beta*Omm+mu*I)^-1;
    V4 = IF2*(beta*Omk+mu*(U-D4));
    %update D
    D1 = D1 - A*U + V1;
    D2 = D2 - U + V2;
    D3 = D3 - U + V3;
    D4 = D4 - U + V4;
    
    if mod(i, 10) == 1
        %object function
        obj = 1/2*norm(A*U - Y, 'fro');
        %primal residual
        res_p = norm([V1; V2; V3; V4] - [A*U; U; U; U], 'fro');
        %dual residual
        res_d = norm([V1; V2; V3; V4] - [V10; V20; V30; V40], 'fro');
        
        if res_p > 10*res_d
            mu = mu*2;D1 = D1/2;D2 = D2/2;D3 = D3/2;D4 = D4/2;
        elseif res_d > 10*res_p
            mu = mu/2;D1 = D1*2;D2 = D2*2;D3 = D3*2;D4 = D4*2;
        end
    end
    sre(i)=20*log10(norm(XT,'fro')/norm(gather(U)-XT,'fro'));
    xreds(i)=norm(gather(U)-XT,'fro');
    obj1(i) = norm(A*U - Y, 'fro');
    %primal residual
    res_p1(i) = norm([V1; V2; V3; V4] - [A*U; U; U; U], 'fro');
    %dual residual
    res_d1(i) = norm([V1; V2; V3; V4] - [V10; V20; V30; V40], 'fro');

    fprintf('i = %d, obj = %.4f,res_p = %.4f, res_d = %.4f, mu = %.1f, SRE=%.2f, ||X-XT||=%.2f, [alpha=%.1e, beta=%.1e]\n',...
        i, obj1(i), res_p1(i), res_d1(i), mu,sre(i),xreds(i),alpha,beta);

    i = i + 1;
end
if i == MaxIter + 1
    display('Maximum iteration reached!');
end
X = gather(U);
end

function Xup=solveV2(X,lambda,imgsize)

%w1=repmat(1./(sum(X.*X)),size(X,1),1);
w1=repmat(1./(sum(X.*X,2)),1,size(X,2));
w2=sw(X,imgsize,'yes');
Xup = soft(X,w1.*w2.*lambda);
end

function W = sw(S,imgsize,filter)
Sg=gather(S);
[p,~]=size(Sg);
nr=imgsize(1);nc=imgsize(2);
if strcmp(filter,'yes')
    Sg=reshape(Sg',nr,nc,p);
    for i=1:p;
        Sg(:,:,i) = wiener2(Sg(:,:,i),[3 3]);
    end
    Sg=reshape(Sg,nr*nc,p)';
    W=1./(abs(Sg)+eps);
else
    W=1./(abs(Sg)+eps);
end
end

function y = soft(x,T)

y = max(abs(x) - T, 0);
y = y./(y+T) .* x;
end
