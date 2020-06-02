clc;
close all
ch0=10;
N0=600;
sorh='h';
keepapp=1;
crit='threshold';
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
   'db41','db42','db43','db44','db45',...        %dbN,N可取1-45； N a positive integer from 1 to 45.
   'coif1','coif2','coif3','coif4','coif5',...   %coiflets 可取N:1-5
   };
s=struct(field,value);
thrx=zeros(1,50);
thr= zeros(1,50);
load X_0.mat;
X=X_0;
level=3;
AMDL=zeros(ch0,50);
jidi=zeros(ch0,1);
N = N0;
lg=log(N)/log(2);

wp=zeros(N,ch0);
%    parpool(40)
%  parfor (ch=1:ch0,40)
 for ch=1:10
     disp(ch)
try
    chanal=ch;
    d=double(X(:,chanal));
    mdl0=zeros(1,50);
    thrx=zeros(1,50);
    thr=zeros(1,50);
    for ji=1:50
        wavename=s(ji).name;
        t0=wpdec(d,level,wavename,crit,1);
        [t,ent]=besttree(t0); 
        wa=read(t,'data');
        wa0=sort(wa,'descend');
        km=length(wa);
        A=zeros(km-1,1);
        for k=1:km-1
            th=wa;
            thr0=wa0(k);
            th(th<thr0)=0;
            M=wa+(-th);
            f1=1.5*k*lg;
            f2=0.5*N*(log(sum(M.*M))/log(2));
            A(k)=f1+f2;
        end
        mdl0(ji)=min(A);
        thrx(ji)=find(A==min(A));
        thr(ji)=wa0(thrx(ji)+1);
    end 
    zu=find( thrx==min( thrx));  
    AMDL(ch,:)=mdl0;
    jidi(ch) =zu(1);
    wavename=s(zu(1)).name;
    th=thr(zu(1));
    [xout,treed,perf0,perfl2] = wpdencmp(d,sorh,level,wavename,crit,th,keepapp);
    wp(:,ch)=xout;
catch
    wp(:,ch)=d;
end
      figure,subplot(2,1,1),plot(d)
      hold on
      subplot(2,1,1),plot(xout)
      legend('raw data','reconstruction data');
      cch0=['wp-' crit];
      cch0=[cch0,'-ch'];
      cch1=[cch0 num2str(ch)];
      cch2=['-base:' wavename];
      nname=[cch1 cch2];
      title(nname);
     
      cha=wp(:,ch)-d;
      subplot(2,1,2),plot(cha);
      title('diff value');
    
    
end
%   save shAMDL AMDL;
%   save shjidi jidi;
%   save shwp wp;
%   h5create('wp-big-sh.h5','/bigwout',[600,200000],'ChunkSize',[600 100000],'Deflate',9);
%   h5write('wp-big-sh.h5','/bigwout',wp);






