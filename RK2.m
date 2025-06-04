function yapprox = RK2(m, x0, xf, y0, h, a2)
% Eulers method: a2=0
% Heuns method: a2=0
% Midpoints method: a2=1
% Ralstons method: a2=2/3
% Ricos method a2=666

% 2nd order Runge-Kutta

x=x0:h:xf;
y=zeros(size(x));
y(1)=y0;

a1=1-a2;

if a2~=0
    q=1/(2*a2),
else
    q=0;
end
for i=1:length(x)-1
    k1=m(x(i),y(i));
    k2=m(x(i)+q*h,y(i)+q*k1*h);
    y(i+1)=y(i)+(a1*k1+a2*k2)*h;
end
yapprox = y;



plot(x,y,'r')
hold on
plot(x,yexact(x),'b')
hold off
xlabel("x")
ylabel("y")
legend("2nd order Runge-Kutta","Exact")

