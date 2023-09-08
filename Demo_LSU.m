clear all
close all
clc


load cuprite_ref
load USGS_1995_Library
A=datalib(:,4:end);
A=A(BANDS,:);

Y=x;
nr=250;nc=191;nb=188;

load TOTAL_X_Cuprite_D5_M7

tic
[X_fast,Xt,Ap,D,iterf,idx,info]=runLSUCuprite(Y,A,30,0.001,0.1,30,[nr,nc],X_lgsu);
time=toc

%%
load USGS_1995_Library.mat
Name=char(names);
Name(idx+3,:)

% fillcolor1=[255, 69, 0]/255; % fillcolors = rand(24, 3);
% fillcolor2=[0, 255, 127]/255;
color1=[221,106,79]/255;%[255,165,16]/255;
color2=[229,168,75]/255;%[12,132,198]/255;
color3=[73,148,196]/255;%[247,77,77]/255;
color4=[77,128,69]/255;%[255 136 132]/255;%[65,183,172 ]/255;
color5=[238,170,154]/255;%[0, 255, 127]/255;
color6=[178,182,182]/255;%[239,239,239]/255;%[153 153 153]/255;
set1=18:23;set2=67:68;set3=81;set4=241:245;set5=288:297;
for i=1:6
    figure
    activeness=info{1,i}{1,2};
    cutoff=info{1,i}{1,5};
    [~,aidx]=sort(activeness,'descend');
    if i==1
        h=plot(set1, activeness(set1),'o','MarkerSize',8,'color',color1) ;
        set(h,'MarkerFaceColor',get(h,'color'));
        hold on
        h=plot(set2, activeness(set2),'s','MarkerSize',8,'color',color2) ;
        set(h,'MarkerFaceColor',get(h,'color'));        
        h=plot(set3, activeness(set3),'d','MarkerSize',8,'color',color3) ;
        set(h,'MarkerFaceColor',get(h,'color'));        
        h=plot(set4, activeness(set4),'<','MarkerSize',8,'color',color4) ;
        set(h,'MarkerFaceColor',get(h,'color'));        
        h=plot(set5, activeness(set5),'>','MarkerSize',8,'color',color5) ;
        set(h,'MarkerFaceColor',get(h,'color'));
        
    else
        idxUpper=info{1,i-1}{1,3};
        [set11,idxA]=intersect(idxUpper,set1);
        if isempty(set11)==0
            h=plot(idxA, activeness(idxA),'o','MarkerSize',8,'color',color1) ;
            set(h,'MarkerFaceColor',get(h,'color'));
            hold on
        end
        [set21,idxA]=intersect(idxUpper,set2);
        if isempty(set21)==0
            h=plot(idxA, activeness(idxA),'s','MarkerSize',8,'color',color2) ;
            set(h,'MarkerFaceColor',get(h,'color'));
        end
        [set31,idxA]=intersect(idxUpper,set3);
        if isempty(set31)==0
            h=plot(idxA, activeness(idxA),'d','MarkerSize',8,'color',color3) ;
            set(h,'MarkerFaceColor',get(h,'color'));
        end
        [set41,idxA]=intersect(idxUpper,set4);
        if isempty(set41)==0
            h=plot(idxA, activeness(idxA),'<','MarkerSize',8,'color',color4) ;
            set(h,'MarkerFaceColor',get(h,'color'));
        end
        [set51,idxA]=intersect(idxUpper,set5);
        if isempty(set51)==0
            h=plot(idxA, activeness(idxA),'>','MarkerSize',8,'color',color5) ;
            set(h,'MarkerFaceColor',get(h,'color'));
        end
    end
    
    hold on
    %xl=1:1:length(activeness);
    h=plot(activeness,'o','MarkerSize',4,'color',color6);
    set(h,'MarkerFaceColor',get(h,'color'));
    plot(repmat(cutoff,1,length(activeness)),'k-','LineWidth',0.5);
    h=xlabel('Dictionary','FontSize', 14);set(h, 'FontName', 'Calibri', 'FontWeight', 'bold')
    h=ylabel('Activeness','FontSize', 14);set(h, 'FontName', 'Calibri', 'FontWeight', 'bold')
    ylim([-50,inf]);xlim([1,length(activeness)])
    if i==6
    h=legend('Alunite','Buddingtonite','Chalcedony','Kaolin/Smect','Montmorillonite','Activeness of libraries','Cutoff value','location','northwest');
    set(h, 'FontName', 'Calibri', 'FontWeight', 'bold','FontSize', 14);
    end
end


%Alunite GDS82 Na82£º              18-23(18-20)--20
%Andradite WS487£º
%Buddingtonite GDS85 D-206£º       67-68--67
%Chalcedony CU91-6A£º              81
%Kaolin/Smect H89-FR-5 30K£º       241-245--244
%Kaolin/Smect KLF508 85%K£º
%Kaolinite KGa-2£º
%Montmorillonite + Illi CM37       288-297--296
%Muscovite IL107£º
%Nontronite NG-1.a£º
%Pyrope WS474£º
%Sphene HS189.3B£º
NameCell=cell(1);
for i=1:498
    NameCell{i,1}=Name(i,:);
end



%% nonnormalized maps
figure
p=5;q=12;
idxMap=[18,67,81,242,296];
nr=250;nc=191;
aa=0.5;bb=0.5;cc=0.7;dd=0.8;ee=0.5;

subplot_tight(p, q, 1,[.003 .003]);
imagesc(reshape(X_sunsal(idxMap(1),:)',nr, nc),[0,aa]);axis image;axis off;colormap(jet)
subplot_tight(p, q, 13,[.003 .003]);
imagesc(reshape(X_sunsal(idxMap(2),:)',nr, nc),[0,bb]);axis image;axis off;colormap(jet)
subplot_tight(p, q, 25,[.003 .003]);
imagesc(reshape(X_sunsal(idxMap(3),:)',nr, nc),[0,cc]);axis image;axis off;colormap(jet)
subplot_tight(p, q, 37,[.003 .003]);
imagesc(reshape(X_sunsal(idxMap(4),:)',nr, nc),[0,dd]);axis image;axis off;colormap(jet)
subplot_tight(p, q, 49,[.003 .003]);
imagesc(reshape(X_sunsal(idxMap(5),:)',nr, nc),[0,ee]);axis image;axis off;colormap(jet)

subplot_tight(p, q, 2,[.003 .003]);
imagesc(reshape(X_sunsaltv(idxMap(1),:)',nr, nc),[0,aa]);axis image;axis off;colormap(jet)
subplot_tight(p, q, 14,[.003 .003]);
imagesc(reshape(X_sunsaltv(idxMap(2),:)',nr, nc),[0,bb]);axis image;axis off;colormap(jet)
subplot_tight(p, q, 26,[.003 .003]);
imagesc(reshape(X_sunsaltv(idxMap(3),:)',nr, nc),[0,cc]);axis image;axis off;colormap(jet)
subplot_tight(p, q, 38,[.003 .003]);
imagesc(reshape(X_sunsaltv(idxMap(4),:)',nr, nc),[0,dd]);axis image;axis off;colormap(jet)
subplot_tight(p, q, 50,[.003 .003]);
imagesc(reshape(X_sunsaltv(idxMap(5),:)',nr, nc),[0,ee]);axis image;axis off;colormap(jet)

subplot_tight(p, q, 3,[.003 .003]);
imagesc(reshape(X_clsunsal(idxMap(1),:)',nr, nc),[0,aa]);axis image;axis off;colormap(jet)
subplot_tight(p, q, 15,[.003 .003]);
imagesc(reshape(X_clsunsal(idxMap(2),:)',nr, nc),[0,bb]);axis image;axis off;colormap(jet)
subplot_tight(p, q, 27,[.003 .003]);
imagesc(reshape(X_clsunsal(idxMap(3),:)',nr, nc),[0,cc]);axis image;axis off;colormap(jet)
subplot_tight(p, q, 39,[.003 .003]);
imagesc(reshape(X_clsunsal(idxMap(4),:)',nr, nc),[0,dd]);axis image;axis off;colormap(jet)
subplot_tight(p, q, 51,[.003 .003]);
imagesc(reshape(X_clsunsal(idxMap(5),:)',nr, nc),[0,ee]);axis image;axis off;colormap(jet)

subplot_tight(p, q, 4,[.003 .003]);
imagesc(reshape(X_drsu(idxMap(1),:)',nr, nc),[0,aa]);axis image;axis off;colormap(jet)
subplot_tight(p, q, 16,[.003 .003]);
imagesc(reshape(X_drsu(idxMap(2),:)',nr, nc),[0,bb]);axis image;axis off;colormap(jet)
subplot_tight(p, q, 28,[.003 .003]);
imagesc(reshape(X_drsu(idxMap(3),:)',nr, nc),[0,cc]);axis image;axis off;colormap(jet)
subplot_tight(p, q, 40,[.003 .003]);
imagesc(reshape(X_drsu(idxMap(4),:)',nr, nc),[0,dd]);axis image;axis off;colormap(jet)
subplot_tight(p, q, 52,[.003 .003]);
imagesc(reshape(X_drsu(idxMap(5),:)',nr, nc),[0,ee]);axis image;axis off;colormap(jet)

subplot_tight(p, q, 5,[.003 .003]);
imagesc(reshape(X_s2wsu(idxMap(1),:)',nr, nc),[0,aa]);axis image;axis off;colormap(jet)
subplot_tight(p, q, 17,[.003 .003]);
imagesc(reshape(X_s2wsu(idxMap(2),:)',nr, nc),[0,bb]);axis image;axis off;colormap(jet)
subplot_tight(p, q, 29,[.003 .003]);
imagesc(reshape(X_s2wsu(idxMap(3),:)',nr, nc),[0,cc]);axis image;axis off;colormap(jet)
subplot_tight(p, q, 41,[.003 .003]);
imagesc(reshape(X_s2wsu(idxMap(4),:)',nr, nc),[0,dd]);axis image;axis off;colormap(jet)
subplot_tight(p, q, 53,[.003 .003]);
imagesc(reshape(X_s2wsu(idxMap(5),:)',nr, nc),[0,ee]);axis image;axis off;colormap(jet)

subplot_tight(p, q, 6,[.003 .003]);
imagesc(reshape(X_qmv(3,:)',nr, nc),[0,aa]);axis image;axis off;colormap(jet)
subplot_tight(p, q, 18,[.003 .003]);
imagesc(reshape(X_qmv(2,:)',nr, nc),[0,bb]);axis image;axis off;colormap(jet)
subplot_tight(p, q, 30,[.003 .003]);
imagesc(reshape(X_qmv(11,:)',nr, nc),[0,cc]);axis image;axis off;colormap(jet)
subplot_tight(p, q, 42,[.003 .003]);
imagesc(reshape(X_qmv(4,:)',nr, nc),[0,dd]);axis image;axis off;colormap(jet)
subplot_tight(p, q, 54,[.003 .003]);
imagesc(reshape(X_qmv(8,:)',nr, nc),[0,ee]);axis image;axis off;colormap(jet)

subplot_tight(p, q, 7,[.003 .003]);
imagesc(reshape(X_mua(idxMap(1),:)',nr, nc),[0,aa]);axis image;axis off;colormap(jet)
subplot_tight(p, q, 19,[.003 .003]);
imagesc(reshape(X_mua(idxMap(2),:)',nr, nc),[0,bb]);axis image;axis off;colormap(jet)
subplot_tight(p, q, 31,[.003 .003]);
imagesc(reshape(X_mua(idxMap(3),:)',nr, nc),[0,cc]);axis image;axis off;colormap(jet)
subplot_tight(p, q, 43,[.003 .003]);
imagesc(reshape(X_mua(idxMap(4),:)',nr, nc),[0,dd]);axis image;axis off;colormap(jet)
subplot_tight(p, q, 55,[.003 .003]);
imagesc(reshape(X_mua(idxMap(5),:)',nr, nc),[0,ee]);axis image;axis off;colormap(jet)

subplot_tight(p, q, 8,[.003 .003]);
imagesc(reshape(X_wcsutv(idxMap(1),:)',nr, nc),[0,aa]);axis image;axis off;colormap(jet)
subplot_tight(p, q, 20,[.003 .003]);
imagesc(reshape(X_wcsutv(idxMap(2),:)',nr, nc),[0,bb]);axis image;axis off;colormap(jet)
subplot_tight(p, q, 32,[.003 .003]);
imagesc(reshape(X_wcsutv(idxMap(3),:)',nr, nc),[0,cc]);axis image;axis off;colormap(jet)
subplot_tight(p, q, 44,[.003 .003]);
imagesc(reshape(X_wcsutv(idxMap(4),:)',nr, nc),[0,dd]);axis image;axis off;colormap(jet)
subplot_tight(p, q, 56,[.003 .003]);
imagesc(reshape(X_wcsutv(idxMap(5),:)',nr, nc),[0,ee]);axis image;axis off;colormap(jet)

load cuprite_SUnCNN_SUresult_1
X3d=reshape(X_pred',nc,nr,size(A,2));
X_suncnn=zeros(nr,nc,size(A,2));
for i=1:nr
    for j=1:nc
        for k=1:size(A,2)
            X_suncnn(i,j,k)=X3d(j,i,k);
        end
    end
end
X_suncnn=reshape(X_suncnn,nr*nc,size(A,2))';
subplot_tight(p, q, 9,[.003 .003]);
imagesc(reshape(X_suncnn(idxMap(1),:)',nr, nc),[0,aa]);axis image;axis off;colormap(jet)
subplot_tight(p, q, 21,[.003 .003]);
imagesc(reshape(X_suncnn(idxMap(2),:)',nr, nc),[0,bb]);axis image;axis off;colormap(jet)
subplot_tight(p, q, 33,[.003 .003]);
imagesc(reshape(X_suncnn(idxMap(3),:)',nr, nc),[0,cc]);axis image;axis off;colormap(jet)
subplot_tight(p, q, 45,[.003 .003]);
imagesc(reshape(X_suncnn(idxMap(4),:)',nr, nc),[0,dd]);axis image;axis off;colormap(jet)
subplot_tight(p, q, 57,[.003 .003]);
imagesc(reshape(X_suncnn(idxMap(5),:)',nr, nc),[0,ee]);axis image;axis off;colormap(jet)

subplot_tight(p, q, 10,[.003 .003]);
imagesc(reshape(X_sslrsu(idxMap(1),:)',nr, nc),[0,aa]);axis image;axis off;colormap(jet)
subplot_tight(p, q, 22,[.003 .003]);
imagesc(reshape(X_sslrsu(idxMap(2),:)',nr, nc),[0,bb]);axis image;axis off;colormap(jet)
subplot_tight(p, q, 34,[.003 .003]);
imagesc(reshape(X_sslrsu(idxMap(3),:)',nr, nc),[0,cc]);axis image;axis off;colormap(jet)
subplot_tight(p, q, 46,[.003 .003]);
imagesc(reshape(X_sslrsu(idxMap(4),:)',nr, nc),[0,dd]);axis image;axis off;colormap(jet)
subplot_tight(p, q, 58,[.003 .003]);
imagesc(reshape(X_sslrsu(idxMap(5),:)',nr, nc),[0,ee]);axis image;axis off;colormap(jet)

subplot_tight(p, q, 11,[.003 .003]);
imagesc(reshape(X_lgsu(idxMap(1),:)',nr, nc),[0,aa]);axis image;axis off;colormap(jet)
subplot_tight(p, q, 23,[.003 .003]);
imagesc(reshape(X_lgsu(idxMap(2),:)',nr, nc),[0,bb]);axis image;axis off;colormap(jet)
subplot_tight(p, q, 35,[.003 .003]);
imagesc(reshape(X_lgsu(idxMap(3),:)',nr, nc),[0,cc]);axis image;axis off;colormap(jet)
subplot_tight(p, q, 47,[.003 .003]);
imagesc(reshape(X_lgsu(idxMap(4),:)',nr, nc),[0,dd]);axis image;axis off;colormap(jet)
subplot_tight(p, q, 59,[.003 .003]);
imagesc(reshape(X_lgsu(idxMap(5),:)',nr, nc),[0,ee]);axis image;axis off;colormap(jet)

subplot_tight(p, q, 12,[.003 .003]);
imagesc(reshape(X_fast(idx==idxMap(1),:)',nr, nc),[0,aa]);axis image;axis off;colormap(jet); %h=colorbar; set(h,'fontsize',14,'FontWeight', 'bold');
subplot_tight(p, q, 24,[.003 .003]);
imagesc(reshape(X_fast(idx==idxMap(2),:)',nr, nc),[0,bb]);axis image;axis off;colormap(jet); %h=colorbar; set(h,'fontsize',14,'FontWeight', 'bold');
subplot_tight(p, q, 36,[.003 .003]);
imagesc(reshape(X_fast(idx==idxMap(3),:)',nr, nc),[0,cc]);axis image;axis off;colormap(jet); %h=colorbar; set(h,'fontsize',14,'FontWeight', 'bold');
subplot_tight(p, q, 48,[.003 .003]);%
imagesc(reshape(X_fast(idx==idxMap(4),:)',nr, nc),[0,dd]);axis image;axis off;colormap(jet); %h=colorbar; set(h,'fontsize',14,'FontWeight', 'bold');
subplot_tight(p, q, 60,[.003 .003]);
imagesc(reshape(X_fast(idx==idxMap(5),:)',nr, nc),[0,ee]);axis image;axis off;colormap(jet); %h=colorbar; set(h,'fontsize',14,'FontWeight', 'bold');


drawnow;
