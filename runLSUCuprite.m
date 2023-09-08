

function [X,Xt,Ap,D,iter,tidx,err]=runLSUCuprite(Y,A,p,alpha,beta,iter,imgsize,XT)


%--------------------------------------------------------------
% Usage:
%
% X = runLSUCuprite(Y, A, p,alpha,beta,iter,imgsize,XT)
%
% ------- Input variables -------------------------------------------
%
%  Y - hyperspectral data matrix with dimensions L(bands) x K(pixels)
%
%  A - spectral library matrix with dimentsions L(bands) x m(spectra)
%
%  parameter.
%            * p - number of active libraies attempted to identify   
%            * alpha - sparsity regularization penalty parameter
%            * beta - ASC constraint penalty parameter
%            * iter - maximum iteration number
%            * imgsize - the image size (row and column) of input bands
%            * xt - the reference gt
%
%  NOTE: PARAMETERS WITH SYMBOL * ARE NECESSARY
%
% ------- Output variables -------------------------------------------
%
% X - the estimated abundance matrix with dimensions p(spectra) x K(pixels)
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




nway=size(A);nway2=size(Y);
D=size(A,2);
tidx=linspace(1,D,D);

% initialize para
opts.maxiter=iter;opts.imgsize=imgsize;opts.verbose=0;
opts.alpha=alpha;opts.beta=beta;opts.xt=XT;
% opts.gamma=0.5;
err=cell(1);
isrun=1;
it=1;



while isrun==1
    Aup=A;
    %muc(it)=shenHyperMutualCoherence(A);
    %% run algorithm
    
    [X,sre] = LSU(Y, A, opts);
    Xf=X;
    %% update dictionary
    Act=sum(X.*X,2);
    [~,aidx]=sort(Act,'descend');
    if D==p
        isrun=0;
    elseif D/2<p
        D=p;
        A(:,aidx(D+1:end))=[];
        tidx(aidx(D+1:end))=[];
        cutoff=Act(aidx(D));
        %update coeff gt
        XT(aidx(D+1:end),:)=[];
        Xf(aidx(D+1:end),:)=[];
    elseif D/2>p
        D=ceil(D/2);
        A(:,aidx(1+D:end))=[];
        tidx(aidx(1+D:end))=[];
        %update coeff gt
        XT(aidx(D+1:end),:)=[];
        Xf(aidx(D+1:end),:)=[];
        cutoff=Act(aidx(D));
    end
    err{it}{1}=X;
    err{it}{2}=Act;
    err{it}{3}=tidx;
    err{it}{4}=sre;
    err{it}{5}=cutoff;
    err{it}{6}=Aup;
    it=it+1;
    %% update para
    opts.beta=min(opts.beta*2,15);
    opts.maxiter=min(opts.maxiter*2,200);
    opts.xt=XT;
    
    %% plot results
    SRE = 20*log10(norm(XT,'fro')/norm(Xf-XT,'fro'));
    
    SPA=length(find(Xf>0.005))/((size(Xf,1)*size(Xf,2)));
    
    RMSE=sqrt(mean2((Xf-XT).^2));
    fprintf('Plot results at current iteration!\n');
    fprintf('SRE=%.2f---Spa=%.2f---RMSE=%.2f \n',SRE,SPA,RMSE);
    
end
Ap=A;
Xt=zeros(nway(2),nway2(2));
Xt(tidx,:)=X;
end
