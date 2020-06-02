clc;
close all
field='name'; 
value={
   'db1','db2','db3','db4','db5',...
   'db6','db7','db8','db9','db10',...
   'db11','db12','db13','db14','db15',...
   'db16','db17','db18','db19','db20',...
   'db21','db22','db23','db24','db25',...
   'db26','db27','db28','db29','db30',...
   'db31','db32','db33','db34','db35',...
   'db36','db37','db38','db39','db40',...
   'db41','db42','db43','db44','db45',...
   'coif1','coif2','coif3','coif4','coif5',...
   
   };
s=struct(field,value);

ch0=1440000;
N0=1029;

load X.mat;
AMDL=zeros(ch0,50);
jidi=zeros(ch0,1);
wavelet=zeros(ch0,N0);
lg=log(N0)/log(2);

parpool(100)
 parfor (ch=1:ch0,100)
% for ch=1:ch0
    disp(ch)
    try
    chanal=ch;  
    d=double(X(:,chanal));
    mdl=zeros(1,50);
    thro=zeros(1,50);    
    for sname=1:50
        wavename=s(sname).name;
%         level=wmaxlev(N0,wavename);
         level=4;
        [wa,~]=wavedec(d,level,wavename);
        wa0=sort(wa); 
        wj= fliplr(wa0);
        mn=zeros(length(wa)-1,1);
        for k=1:(length(wa)-1)
            th=wa;
            thr0=wj(length(wj)+1-k);
            th(th<thr0)=0;
            M=wa-th;
            f1=1.5*k*lg;
            f2=0.5*N0*(log(sum(M.*M))/log(2));
            mn(k)=f1+f2;%
        end
        thrx=find(mn==min(mn));
        mdl(sname)=min(mn);
        thro(sname)=wj(length(wj)-thrx+1);
    end
    we=find(mdl(:)==min(mdl(:)));
    jidi(ch,1)=we;
    bestwave=s(we).name;
    wavename=bestwave;          
    [wa,L]=wavedec(d,level,wavename);
    thrx=find(mn==min(mn));
    thr=wj(length(wj)-thrx);
    wa(wa<thr)=0;
    wavelet(ch,:)=fix(waverec(wa,L,wavename));
    AMDL(ch,:)=mdl;
    catch
        wavelet(ch,:)=d;
    end
    
%     figure,subplot(2,1,1),plot(d)
%     hold on
%     subplot(2,1,1),plot( wavelet(ch,:))
%     legend('raw data','reconstruction data');
%     cch1=['w-ch' num2str(ch)];
%     cch2=['-base:' wavename];
%     nname=[cch1 cch2];
%     title(nname);
%     
%     cha=wavelet(ch,:)-d';
%     subplot(2,1,2),plot(cha);
%     title('diff value');
%      saveas(gcf,[cch1 '.jpg']);
end
save wlAMDL AMDL;
save wljidi jidi;
h5create('w-lownoise.h5','/bigwout',[ch0,N0],'Datatype','int16','ChunkSize',[10000,N0],'Deflate',9);
h5write('w-lownoise.h5','/bigwout',wavelet);




