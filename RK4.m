function yapprox = RK4(m, x0, xf, y0, h, a2)

x=x0:h:xf;
y=zeros(size(x));
y(1)=y0;

for i=1:length(x)-1
    k1=m(x(i),y(i));
    k2=m(x(i)+h/2,y(i)+k1*h/2);
    k3=m(x(i)+h/2,y(i)+k2*h/2);
    k4=m(x(i)+h,y(i)+k3*h);
    y(i+1)=y(i)+((1/6*k1+(1/2)*k2+(1/2)*k3+(1/6)*k4))*h;
end
yapprox = y;

