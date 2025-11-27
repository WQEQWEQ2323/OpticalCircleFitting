function y = struve(dimension,x)
% y=y+a abs(a)先增大，后减小
%x的最大值不能过大，否则可能a还来不及收敛就超出了数值上限
j=1e3;
a=zeros(1,length(x));
y=zeros(1,length(x));
a=1/gamma(3/2)/gamma(dimension+3/2)*(x/2).^(dimension+1);
for i=2:j
    if max(abs(a))<1e-3
        break;
    end
    y=y+a;
    a=a.*(x.^2)/(i-1/2)/(dimension+i-1/2)*(-1/4);
end