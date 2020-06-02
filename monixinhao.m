t=1:1029;
x=exp(-0.6*t);
x=-30*x+967;
x1=x;
for i=1:1029
    X(i)=x1(1029-i+1);
end
x2=exp(-0.3*t);
x2=-30*x2+967;
x3=x2;
for i=1:200
    XX(i)=X(829+i);
end
for i=201:1029
    XX(i)=x3(i-200);
end

X2=XX;
y=wgn(1,1029,2);
figure,plot(y);
X=X2;
diej=X+y;
figure,plot(diej);
 as=SNR_singlech(diej,XX);


