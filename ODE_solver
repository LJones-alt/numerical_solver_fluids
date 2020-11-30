function y = ode(eta,f) 
% f()  f = f(1)  f' = f(2) f" = f(3) 
b=0.05; %commnent out for section 1
y = zeros(3,1); %section 1 and 2
y = zeros(5,1); %comment out unless compressible 3 only
% y() is array of the derivatives of f() 
y(1) = f(2);           
y(2) = f(3); 
y(3) = (-f(1)*f(3)); %no transpiration
y(3) = (-f(1)*f(3))-b*(1-(f(2)*f(2))); %FK equation with beta values- comment out as apropriate
y(3) = (-f(1)*f(3))-b*(1-(f(2)*f(2)))-f(5)); %for use in compressible flow Sw values, comment out as appropriate
y(4) = f(4); %comment out for blasis and FK
y(5) = -f(1)*f(4); %commmet out for blasius and FK
