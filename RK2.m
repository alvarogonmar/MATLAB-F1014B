function yapprox = RK2()

% 2nd order Runge-Kutta
y0=10
x0=0
xf=3
h=1

x=x0:h:xf
y=zeros(size(x))
y(1)=y0
err=zeros(size(x))

a2=1
a1=1-a2
q=1/(2*a2)


for i=1:length(x)-1
    k1=m(x(i),y(i));
    k2=m(x(i)+q*h,y(i)+q*k1*h);
    y(i+1)=y(i)+(a1*k1+a2*k2)*h;

    err(i+1)=abs((yexact(x(i+1))-y(i+1))/yexact(x(i+1)))*100;
end


plot(x,y,'r')
hold on
plot(x,yexact(x),'b')
hold off
xlabel("x")
ylabel("y")
legend("2nd order Runge-Kutta","Exact")

