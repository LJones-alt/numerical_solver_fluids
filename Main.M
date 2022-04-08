%%
%ODE solver for Blausius equns (05/12/19 . V1.5 Laura Jones date: 20/12/19
%%
guess_1 = 0.1;
guess_2 = 0.95; %make sufficiently large to capture above answer
guess_3 = -1; %these get removed for section 1 and 2
guess_4 = 0.5;
tol = 0.001; %tolerance of the answer- set as desired
error = 1; %`set guess1 error 
s_error = 1;  %set to zero for incompressible
count = 1; %count iterations
 Sw =1;  %comment out for section 1 and 2
 
v_w = 0; %values of Vw set to 0 for no-transpiration
 
%%
%conditional while loop
     while error > tol && s_error >tol  %while error is bigger than specified tolerance, do this      
        [eta,y]=ode45(@ode,[0 10],[v_w;0; guess_1]); %guess_3 ; Sw]); %use ode45 to solve system using ode function created
        [eta_2,y_2]=ode45(@ode,[0 10],[v_w; 0;guess_2]); % guess_4; Sw]); %same as above
        tau_1 = y(end,2);  %extract last value eta= inf 
        tau_2 = y_2(end,2); %same as above with guess 2       
        S_1 = y(end,5);
        S_2 = y_2(end,5);
        error_1 = abs(1-tau_1);  %find error from known value of f''(inf) = 0
        error_2 = abs(1-tau_2); %same as above for second guess
        error_3= abs(S_1); %comment out for section 1 and 2
        error_4 = abs(S_2);  %comment out for section 1 and 2    
        midpoint = 0.5*(guess_1 + guess_2); %bisection to find next guess
        midpoint_s = 0.5*(guess_3 + guess_4);%comment out for section 1 and 2
        if error_1 > error_2  %if guess 2 was closer, make this into new guess 1
            tau_1 = tau_2;  %make better guess the one to use 
            tau= y_2(:,2);
            error = error_2; %transfer error over
            guess_1= guess_2;
            guess_2 = midpoint;
        elseif error_1< error_2  %if guess 1 was better, then the error for the loop is defined by guess 1
            error = error_1;
            tau=y(:,2);
            guess_2 = midpoint;
        end      
        %include this loop for section 3, otherwise comment out
        if error_3 > error_4
            S_1 = S_2;
            s_error = error_4;
            guess_3 = guess_4;
            guess_4 = midpoint_s;
           
        elseif error_3 < error_4
            s_error = error_3;
            guess_4 = midpoint_s;
        end
         disp("current error  " + error + "  " + s_error); %comment out if error convergence not needed
         count = count +1;
         disp("Current iteration   " + count);   %displays current iteration number
     end%end of while loop
    tau_b = y(1,3);
 
  %%
  %plot results as the values change
    hold on
     plot(eta_2,y_2) %use one with smaller error
     plot(eta,y)
 
  title("Solution with B=0.05 and range of Vw ") %change as according
 xlabel("eta")
 ylabel("Value of Function f,f',f''") %include S and S' in 3rd section
 legend('f','f-dash','F-dash-dash'); %change as according
 
%%
%post processing
%set up of variables
x=zeros(1001,11);
eta=zeros(1001,11);
U_inf = 1;
L = 10;
mu = 1.789E-5;
rho = 1.225;
nu = mu/rho;
A = sqrt(nu/U_inf);
h = 0.01;
 
%delta- finds bounday layer thickness
for i = 1:length(x)
delta(i) = 5*sqrt(x(i))*A;
end
plot(x, delta, 'color', 'black');
position = [1 5 10];
for j = 1:length(position)
y(j) = eta*sqrt(position(j))*A;
plot(y2+position(j), y{j})
end
%drag coeff - calculate Cf over plate cumulativly
y_3=y(:,3);
for i = 1 : length(x)
    tau_x(i) = mu*U_inf*y_3(1)*sqrt(U_inf/2/nu/x(i));
end
for i = 1: length(x)
cfx(i) = tau_x(i)*2/rho/U_inf^2;
% cfxt(i) = 0.664/sqrt(rho*U_inf*x(i)/mu); %incompressible
end
 
%skin friction- calculate skin friction over plate
CF_T = 0;
dx = x(2)-x(1);
for i = 2 : length(x)
CF_T = CF_T + 1/L*cfx(i)*dx;
end
CDF_T = CF_T*2;
CF = 0;
for i = 2 : length(x)
CF(i) = 1/L*cfx(i)*dx;
CF(i) = CF(i-1)+CF(i);
CDF(i) = CF(i)*2;
D(i) = 1/2*rho*U_inf^2*L*CDF(i);
end
