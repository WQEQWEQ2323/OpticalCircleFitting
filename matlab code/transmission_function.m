function result=transmission_function(r,R)
%子函数序号#7
%result=zeros(1,length(r));
result=pi*R/2*(besselj(1,r*R).*StruveH0(r*R)-besselj(0,r*R).*StruveH1(r*R))./r;
result(r==0)=0;
end