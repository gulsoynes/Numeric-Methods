%% Definiton
% This code contains algorithms for numerical methods. It is written as
% part of Numerical Methods course by following the book of

% Chapra,S.C.,and Canale,R.P., Numerical Methods for Engineers, 7th ed.,...
% McGraw-Hill Education, 2015.

% f = function
% Step = structure for solution. It contains percentage error.
% err_limit = desired error% level

%% CLOSED METHODS
%% 1.GRAPHICAL APPROACH
clear
clc
syms x
%f=5*x^3-5*x^2+6*x-2;
f=log(x^2)-0.7;
Desired=0;
g(x)=f-Desired;
fplot(g(x))
grid on
eq=g(x)==0;
[sol]=vpasolve(eq,x);

%% 2.BISECTION METHOD
clear
clc
syms x
f=x^3-6*x^2+11*x-6.1;
Desired=0;
err_limit=1;
g(x)=f-Desired;

eq=g(x)==0;

[sol]=vpasolve(eq,x);
x_r_true=(double([sol]));
Step(1).x_r_t=x_r_true;
x_upper=3.5;
x_lower=2.5;
Step(2).x_r_t=x_r_true;
xr=(x_upper+x_lower)/2;
for i=1:20
    
    Step(i+1).xup=x_upper;
    Step(i+1).xlow=x_lower;
    
    a=double(g(x_upper));
    Step(i+1).g_up=a;
    b=double(g(x_lower));
    Step(i+1).g_lw=b;
    
    xr=(x_upper+x_lower)/2;
    Step(i+1).x_r=xr;
    c=double(g(xr));
    
    Step(i+1).g_x_r=c;

    if g(xr)*g(x_lower)<0
        x_upper=xr;
    else
        if g(xr)*g(x_lower)>0
            x_lower=xr;
            
        else
            if g(xr)*g(x_lower)==0
                xr=xr;
            end
            
        end
    end
    
    error=(Step(i+1).x_r-Step(i).x_r)/Step(i+1).x_r*100;
    Step(i+1).err=abs(error);
    error_t=(x_r_true-Step(i+1).x_r)/x_r_true*100;
    Step(i+1).trerror=abs(error_t);
    
    if Step(i+1).err<=err_limit
        break
    end
end

%% FALSE POSITION METHOD
clear
clc
syms x
%f=(667.38/x)*(1-exp(-0.146843*x));
%f=x^10-1;
%f=log(x^2)-0.7;
f=x*exp((-2*x))+x^2+2;
Desired=0;
err_limit=0.1;
g(x)=f-Desired;
eq=g(x)==0;
[sol]=vpasolve(eq,x);
x_r_true=(double([sol]));
Step(1).x_r_t=x_r_true;
x_upper=0;
x_lower=-1;
a=double(g(x_upper));
b=double(g(x_lower));
x_r_true=x_r_true(x_lower<real(x_r_true)&real(x_r_true)<x_upper&imag(x_r_true)==0);
Step(2).x_r_t=x_r_true;
for i=1:20
    Step(i+1).xup=x_upper;
    Step(i+1).xlow=x_lower;
    a=double(g(x_upper));
    Step(i+1).g_up=a;
    b=double(g(x_lower));
    Step(i+1).g_lw=b;
    xr=x_upper-(a*(x_lower-x_upper))/(b-a);
    Step(i+1).x_r=xr;
    c=double(g(xr));
    Step(i+1).g_x_r=c;
    %if c==0
    %   break;
    %end
    if g(xr)*g(x_lower)<0
        x_upper=xr;
    else
        if g(xr)*g(x_lower)>0
            x_lower=xr;
            
        else
            if g(xr)*g(x_lower)==0
                xr=xr;
            end
        end
    end
    error=(Step(i+1).x_r-Step(i).x_r)/Step(i+1).x_r;
    Step(i+1).err=abs(error);
    error_t=(x_r_true-Step(i+1).x_r)/x_r_true*100;
    Step(i+1).trerror=abs(error_t);
    if Step(i+1).err<err_limit
        break
    end
end

%% OPEN METHODS
%% 1.SIMPLE FIXED ITERATION
clear
clc
syms x
f(x)=x-4*cos(1/sqrt(x));
g(x)=4*cos(1/sqrt(x));
dif=diff(g(x));
x_o=2.5;
err_limit=0.1;
eq=f(x)==0;
[sol]=vpasolve(eq,x);
x_r_true=(double([sol]));
Step(1).x_r=x_o;
dif(x)=diff(g(x));
i=1;
while 1
    if abs(double(dif(Step(i).x_r)))<1  %converge
        a=double(g(Step(i).x_r));
        Step(i+1).x_r=a;
        error=(Step(i+1).x_r-Step(i).x_r)/Step(i+1).x_r*100;
        Step(i+1).err=abs(error);
        error_t=(x_r_true-Step(i+1).x_r)/x_r_true*100;
        if Step(i+1).err<err_limit
            break
        end       
    else  %diverge
        break
    end
    i=i+1;
end
%% 2.NEWTON-RAPSHSON METHOD
clear
clc
syms x
err_limit=1e-16;
f(x)=x^2-5;
x_o=2;
eq=f(x)==0;
[sol]=vpasolve(eq,x);
x_r_true=(double([sol]))
Step(1).x_r_t=x_r_true; %True root of function
dif(x)=diff(f(x))
Step(1).x_r=x_o;
for i=1:100
        a=double(f(Step(i).x_r));
        Step(i).f=a;
        b=double(dif(Step(i).x_r));
        Step(i).df=b;
        Step(i+1).x_r=Step(i).x_r-(a/b);
        error=(Step(i+1).x_r-Step(i).x_r)/Step(i+1).x_r*100;
        Step(i+1).err=abs(error);
        error_t=(x_r_true-Step(i+1).x_r)/x_r_true*100;
        Step(i+1).trerror=abs(error_t);
        if Step(i+1).err<err_limit
            break
        end
end
%% 3.SECANT METHOD
clear
clc
syms x
err_limit=10^-6;
f(x)=exp(-x)-x;
x_o1=0;  %xi-1
x_o=1;   %xi
eq=f(x)==0;
[sol]=vpasolve(eq,x);
x_r_true=(double([sol]));
Step(1).x_r_t=x_r_true;
Step(2).x_r_t=x_r_true;
Step(1).x_r=x_o1;  %xi-1
Step(2).x_r=x_o;   %xi
for i=1:20
        a=double(f(Step(i).x_r));
        Step(i).f=a; %xi-1
        b=double(f(Step(i+1).x_r));
        %Step(i+1).fo=b;  %xi
        Step(i+2).x_r=Step(i+1).x_r-((b*(Step(i).x_r-Step(i+1).x_r)/(a-b)));
        error=(Step(i+2).x_r-Step(i+1).x_r)/Step(i+2).x_r;
        Step(i+2).err=abs(error);
        error_t=(x_r_true-Step(i+2).x_r)/x_r_true*100;
        Step(i+2).trerror=abs(error_t);
        if Step(i+2).err<err_limit
            break
        end       
end

%% 4.MODIFIDED SECANT METHOD
clear
clc
syms x
err_limit=10^-6;
%f(x)=x*exp(-2*x)+x^2+2;
f(x)=exp(-x)-x;
x_o=1;  %xi-1
delx=0.01;
eq=f(x)==0;
[sol]=vpasolve(eq,x);
x_r_true=(double([sol]));
Step(1).x_r_t=x_r_true;
Step(2).x_r_t=x_r_true;
Step(1).x_r=x_o;  %xi-1
for i=1:20
        a=double(f(Step(i).x_r));
        Step(i).f=a; %xi-1
        b=double(f(Step(i).x_r+delx));
        Step(i).fo=b;  %xi
        Step(i+1).x_r=Step(i).x_r-((a*(Step(i).x_r*delx)/(b-a)));
        error=(Step(i+1).x_r-Step(i).x_r)/Step(i+1).x_r*100;
        Step(i+1).err=abs(error);
        error_t=(x_r_true-Step(i+1).x_r)/x_r_true*100;
        Step(i+1).trerror=abs(error_t);
        if Step(i+1).err<err_limit
            break
        end       
end