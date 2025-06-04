function yapprox = RK4(m, x0, xf, y0, h, a2)

x=x0:h:xf;
y=zeros(size(x));
y(1)=y0;

for i=1:length(x)-1
    k1=m(x(i),y(i));
    k2=m(x(i)+q*h,y(i)+q*k1*h);
    y(i+1)=y(i)+(a1*k1+a2*k2)*h;
end
yapprox = y;

