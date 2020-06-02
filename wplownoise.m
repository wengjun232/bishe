clc;
close all
ch0=1440000;
N0=1029;
sorh='h';
keepapp=1;
crit='shannon';
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
load X.mat;

level=4;
AMDL=zeros(ch0,50);
jidi=zeros(ch0,1);
N = N0;
lg=log(N)/log(2);

wp=zeros(N,ch0);
parpool(100)
parfor (ch=1:ch0,100)
    disp(ch)
try
    chanal=ch;
    d=double(X(:,chanal));
    mdl0=zeros(1,50);
    thrx=zeros(1,50);
    thr=zeros(1,50);
    for ji=1:50
        wavename=s(ji).name;
        t0=wpdec(d,level,wavename,crit);
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
        thr(ji)=wa0(thrx(ji));
    end 
    zu=find( thrx==min( thrx));  
    AMDL(ch,:)=mdl0;
    jidi(ch) =zu(1);
    wavename=s(zu(1)).name;
    th=thr(zu(1));
    [xout,treed,perf0,perfl2] = wpdencmp(d,sorh,level,wavename,crit,th,keepapp);
    wp(:,ch)=fix(xout);
%     huitu(d,xout,ch,wavename);
catch
    wp(:,ch)=fix(d);
end
end
 save shAMDL AMDL;
  save shjidi jidi;
  
  h5create('wp-low-sh1.h5','/bigwout',[N0,ch0],'ChunkSize',[N0 10000],'Deflate',9);
  h5write('wp-low-sh1.h5','/bigwout',wp);





