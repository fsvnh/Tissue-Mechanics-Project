%%%PROGRAM CREATED BY FAH SYSAVANH

%READING FILE NAMES
%data contains load and extension
aorta2data = xlsread('aorta2tensiledata.xlsx');
hydrogel5data = xlsread('hydrogel5v1tensiledata.xlsx');
hydrogel10data = xlsread('hydrogel10v1tensiledata.xlsx');
hydrogel15data = xlsread('hydrogel15v1tensiledata.xlsx');

lambdaaortadata = aorta2data(:,7);
hydrogel5lambdadata = hydrogel5data(:,8);
hydrogel10lambdadata = hydrogel10data(:,8);
hydrogel15lambdadata = hydrogel15data(:,8);

finlambdaaorta = (30.68+lambdaaortadata(10274))/30.68;
finhydrogel5lambda = (20.59+hydrogel5lambdadata(1445))/20.59;
finhydrogel10lambda = (19.11+hydrogel10lambdadata(1157))/19.11;
finhydrogel15lambda = (21.75+hydrogel15lambdadata(1994))/21.75;

%FINDING DEFORMATION GRADIENT

F_aorta = [1/sqrt(finlambdaaorta), 0, 0;
            0, finlambdaaorta, 0;
            0, 0, 1/sqrt(finlambdaaorta)];
        

       
F_h5 = [1/sqrt(finhydrogel5lambda), 0, 0;
            0, finhydrogel5lambda, 0;
            0, 0, 1/sqrt(finhydrogel5lambda)];

 
 F_h10 = [1/sqrt(finhydrogel10lambda), 0, 0;
            0, finhydrogel10lambda, 0;
            0, 0, 1/sqrt(finhydrogel10lambda)];
        
 

 F_h15 = [1/sqrt(finhydrogel15lambda), 0, 0;
            0, finhydrogel15lambda, 0;
            0, 0, 1/sqrt(finhydrogel15lambda)];
        
      
 
 disp(F_aorta);
 disp(F_h5);
 disp(F_h10); 
 disp(F_h15);   
 
%FINDING RIGHT CAUCHY-GREEN TENSOR (C = F_T*F)
F_aorta_T = transpose(F_aorta);
F_h5_T = transpose(F_h5);
F_h10_T = transpose(F_h10);
F_h15_T = transpose(F_h15);

C_aorta = F_aorta_T*F_aorta;
C_h5 = F_h5_T*F_h5;
C_h10 = F_h10_T*F_h10;
C_h15 = F_h15_T*F_h15;

disp(C_aorta);
disp(C_h5);
disp(C_h10);
disp(C_h15);

%FINDING GREEN STRAIN: E = 0.5*(C-I)
I_matrix = [1, 0, 0;
            0, 1, 0;
            0, 0, 1];
        
E_aorta = 0.5*(C_aorta - I_matrix);
E_h5 = 0.5*(C_h5 - I_matrix);
E_h10 = 0.5*(C_h10 - I_matrix);
E_h15 = 0.5*(C_h15 - I_matrix);

disp(E_aorta)
disp(E_h5)
disp(E_h10)
disp(E_h15)

%FINDING INFINITESIMAL STRAIN (ep = 0.5*(F+F_t) - I

ep_aorta = 0.5*(F_aorta + F_aorta_T) - I_matrix;
ep_h5 = 0.5*(F_h5 + F_h5_T) - I_matrix;
ep_h10 = 0.5*(F_h10 + F_h10_T) - I_matrix;
ep_h15 =0.5*(F_h15 + F_h15_T) - I_matrix;

disp(ep_aorta)
disp(ep_h5)
disp(ep_h10)
disp(ep_h15)

%STRAIN ENERGY MODEL

%AORTA 2 TEST
syms c c1 c2 c3 E22 E11
E11test = aorta2data(:,10);
E22test = aorta2data(:,11);
S11test = aorta2data(:,12);
S22test = aorta2data(:,13);

plot(E22test(1:4592),S22test(1:4592))
plot(E11test(1:4592),S11test(1:4592))

n = 4592;

obj_fun = @(x)...
+ sum(((x(1).*exp(E11test(1:n).^2.*x(2)+E22test(1:n).^2.*x(3)+E11test(1:n).*E22test(1:n).*x(4).*2.0).*(E11test(1:n).*x(2).*2.0+E22test(1:n).*x(4).*2.0).*(1.0./2.0))-S11test(1:n)).^2+((x(1).*exp(E11test(1:n).^2.*x(2)+E22test(1:n).^2.*x(3)+E11test(1:n).*E22test(1:n).*x(4).*2.0).*(E11test(1:n).*x(4).*2.0+E22test(1:n).*x(3).*2.0).*(1.0./2.0))-S22test(1:n)).^2);

% Minimization of objective function using fmincon
x0=[1,4,1,1]; % initial guesses
lb=[0,-inf,-inf,-inf]; % lower bound for parameters 
ub=[inf,inf,inf,inf]; % upper bound for parameters
A=[0,-1,0,1;0,0,-1,1]; % Variables to enforce inequality constraint c1>c3 
b=[0;0]; % Variables to enforce inequality constraint c2>c3
options = optimoptions('fmincon','Display','iter','Algorithm','sqp'); %Options to control execution of fmincon
% options = optimoptions('fmincon','Display','iter','Algorithm','sqp','MaxIterations',10000,'Ma xFunEval',100000,'StepTolerance',1e-20);
x = fmincon(obj_fun,x0,A,b,[],[],lb,ub,[],options);

c = x(1);
c1 = x(2);
c2 = x(3);
c3 = x(4);

Q = (c1*E11test.*E11test)+(c2*E22test.*E22test)+(2*c3*E11test.*E22test);   %EXP TERM FOR DISEASED STATE

%EQUATION FOR S11 AND S22 IN DISEASED STATE STATE
S11_a2_eqn = (c/2)*((2*c1*E11test)+(2*c3*E22test)).*exp(Q);
S22_a2_eqn = (c/2)*((2*c3*E11test)+(2*c2*E22test)).*exp(Q);

E11points = linspace(0,min(E11test)*-1,4592);
E22points = linspace(0,max(E22test),4592);

figure()
hold on
scatter(E11test(1:4592),S11test(1:4592))
scatter(E22test(1:4592),S22test(1:4592))
plot(E11points,S11_a2_eqn(1:4592),E22points,S22_a2_eqn(1:4592))
hold off
title('Aorta');
xlabel('Green Strain (E) Points');
ylabel('2nd PK Stress (Pa)');
legend('S11 Experimental Values','S22 Experimental Values','S11 Modeled Equation','S22 Modeled Equation');

%HYDROGEL 5% CONCENTRACTION

%FROM HYDROGEL DATA FILE
E11h5data = hydrogel5data(:,11);
E22h5data = hydrogel5data(:,12);
S11h5data = hydrogel5data(:,13);
S22h5data = hydrogel5data(:,14);

plot(E22h5data,S22h5data)

n = 1445;

obj_fun = @(x)...
+ sum(((x(1).*exp(E11h5data(1:n).^2.*x(2)+E22h5data(1:n).^2.*x(3)+E11h5data(1:n).*E22h5data(1:n).*x(4).*2.0).*(E11h5data(1:n).*x(2).*2.0+E22h5data(1:n).*x(4).*2.0).*(1.0./2.0))-S11h5data(1:n)).^2+((x(1).*exp(E11h5data(1:n).^2.*x(2)+E22h5data(1:n).^2.*x(3)+E11h5data(1:n).*E22h5data(1:n).*x(4).*2.0).*(E11h5data(1:n).*x(4).*2.0+E22h5data(1:n).*x(3).*2.0).*(1.0./2.0))-S22h5data(1:n)).^2);

% Minimization of objective function using fmincon
x0=[1,4,1,1]; % initial guesses
lb=[0,-inf,-inf,-inf]; % lower bound for parameters 
ub=[inf,inf,inf,inf]; % upper bound for parameters
A=[0,-1,0,1;0,0,-1,1]; % Variables to enforce inequality constraint c1>c3 
b=[0;0]; % Variables to enforce inequality constraint c2>c3
options = optimoptions('fmincon','Display','iter','Algorithm','sqp'); %Options to control execution of fmincon
% options = optimoptions('fmincon','Display','iter','Algorithm','sqp','MaxIterations',10000,'Ma xFunEval',100000,'StepTolerance',1e-20);
x = fmincon(obj_fun,x0,A,b,[],[],lb,ub,[],options);

c = x(1);
c1 = x(2);
c2 = x(3);
c3 = x(4);

Q = (c1*E11h5data.*E11h5data)+(c2*E22h5data.*E22h5data)+(2*c3*E11h5data.*E22h5data);   %EXP TERM FOR DISEASED STATE

%EQUATION FOR S11 AND S22 IN DISEASED STATE STATE
S11h5eqn = (c/2)*((2*c1*E11h5data)+(2*c3*E22h5data)).*exp(Q);
S22h5eqn = (c/2)*((2*c3*E11h5data)+(2*c2*E22h5data)).*exp(Q);

E11h5points = linspace(0,min(E11h5data)*-1,1445);
E22h5points = linspace(0,max(E22h5data),1445);

figure()
hold on
scatter(E11h5data,S11h5data);
scatter(E22h5data,S22h5data);
plot(E11h5points,S11h5eqn,E22h5points,S22h5eqn);
hold off
title('5% Gelatin Hydrogel');
xlabel('Green Strain (E) Points');
ylabel('2nd PK Stress (Pa)');
legend('S11 Experimental Values','S22 Experimental Values','S11 Modeled Equation','S22 Modeled Equation');

%HYDROGEL 10%
syms c c1 c2 c3 E22 E11

E11h10data = hydrogel10data(:,11);
E22h10data = hydrogel10data(:,12);
S11h10data = hydrogel10data(:,13);
S22h10data = hydrogel10data(:,14);

n = 1157;

obj_fun = @(x)...
+ sum(((x(1).*exp(E11h10data(1:n).^2.*x(2)+E22h10data(1:n).^2.*x(3)+E11h10data(1:n).*E22h10data(1:n).*x(4).*2.0).*(E11h10data(1:n).*x(2).*2.0+E22h10data(1:n).*x(4).*2.0).*(1.0./2.0))-S11h10data(1:n)).^2+((x(1).*exp(E11h10data(1:n).^2.*x(2)+E22h10data(1:n).^2.*x(3)+E11h10data(1:n).*E22h10data(1:n).*x(4).*2.0).*(E11h10data(1:n).*x(4).*2.0+E22h10data(1:n).*x(3).*2.0).*(1.0./2.0))-S22h10data(1:n)).^2);

% Minimization of objective function using fmincon
x0=[1,4,1,1]; % initial guesses
lb=[0,-inf,-inf,-inf]; % lower bound for parameters 
ub=[inf,inf,inf,inf]; % upper bound for parameters
A=[0,-1,0,1;0,0,-1,1]; % Variables to enforce inequality constraint c1>c3 
b=[0;0]; % Variables to enforce inequality constraint c2>c3
options = optimoptions('fmincon','Display','iter','Algorithm','sqp'); %Options to control execution of fmincon
% options = optimoptions('fmincon','Display','iter','Algorithm','sqp','MaxIterations',10000,'Ma xFunEval',100000,'StepTolerance',1e-20);
x = fmincon(obj_fun,x0,A,b,[],[],lb,ub,[],options);

c = x(1);
c1 = x(2);
c2 = x(3);
c3 = x(4);

Q = (c1*E11h10data.*E11h10data)+(c2*E22h10data.*E22h10data)+(2*c3*E11h10data.*E22h10data);   %EXP TERM FOR DISEASED STATE

%EQUATION FOR S11 AND S22 IN DISEASED STATE STATE
S11h10eqn = (c/2)*((2*c1*E11h10data)+(2*c3*E22h10data)).*exp(Q);
S22h10eqn = (c/2)*((2*c3*E11h10data)+(2*c2*E22h10data)).*exp(Q);

E11h10points = linspace(0,min(E11h10data)*-1,1157);
E22h10points = linspace(0,max(E22h10data),1157);

figure()
hold on
scatter(E11h10data,S11h10data);
scatter(E22h10data,S22h10data);
plot(E11h10points,S11h10eqn,E22h10points,S22h10eqn);
hold off
title('10% Gelatin Hydrogel');
xlabel('Green Strain (E) Points');
ylabel('2nd PK Stress (Pa)');
legend('S11 Experimental Values','S22 Experimental Values','S11 Modeled Equation','S22 Modeled Equation');

%HYDROGEL 15%
E11h15data = hydrogel15data(:,11);
E22h15data = hydrogel15data(:,12);
S11h15data = hydrogel15data(:,13);
S22h15data = hydrogel15data(:,14);

n = 1500;

obj_fun = @(x)...
+ sum(((x(1).*exp(E11h15data(1:n).^2.*x(2)+E22h15data(1:n).^2.*x(3)+E11h15data(1:n).*E22h15data(1:n).*x(4).*2.0).*(E11h15data(1:n).*x(2).*2.0+E22h15data(1:n).*x(4).*2.0).*(1.0./2.0))-S11h15data(1:n)).^2+((x(1).*exp(E11h15data(1:n).^2.*x(2)+E22h15data(1:n).^2.*x(3)+E11h15data(1:n).*E22h15data(1:n).*x(4).*2.0).*(E11h15data(1:n).*x(4).*2.0+E22h15data(1:n).*x(3).*2.0).*(1.0./2.0))-S22h15data(1:n)).^2);

% Minimization of objective function using fmincon
x0=[1,4,1,1]; % initial guesses
lb=[0,-inf,-inf,-inf]; % lower bound for parameters 
ub=[inf,inf,inf,inf]; % upper bound for parameters
A=[0,-1,0,1;0,0,-1,1]; % Variables to enforce inequality constraint c1>c3 
b=[0;0]; % Variables to enforce inequality constraint c2>c3
options = optimoptions('fmincon','Display','iter','Algorithm','sqp'); %Options to control execution of fmincon
% options = optimoptions('fmincon','Display','iter','Algorithm','sqp','MaxIterations',10000,'Ma xFunEval',100000,'StepTolerance',1e-20);
x = fmincon(obj_fun,x0,A,b,[],[],lb,ub,[],options);

c = x(1);
c1 = x(2);
c2 = x(3);
c3 = x(4);

Q = (c1*E11h15data.*E11h15data)+(c2*E22h15data.*E22h15data)+(2*c3*E11h15data.*E22h15data);   %EXP TERM FOR DISEASED STATE

%EQUATION FOR S11 AND S22 IN DISEASED STATE STATE
S11h15eqn = (c/2)*((2*c1*E11h15data)+(2*c3*E22h15data)).*exp(Q);
S22h15eqn = (c/2)*((2*c3*E11h15data)+(2*c2*E22h15data)).*exp(Q);

E11h15points = linspace(0,min(E11h15data)*-1,1994);
E22h15points = linspace(0,max(E22h15data),1994);

figure()
hold on
scatter(E11h15data,S11h15data);
scatter(E22h15data,S22h15data);
plot(E11h15points,S11h15eqn,E22h15points,S22h15eqn);
hold off
title('15% Gelatin Hydrogel');
xlabel('Green Strain (E) Points');
ylabel('2nd PK Stress (Pa)');
legend('S11 Experimental Values','S22 Experimental Values','S11 Modeled Equation','S22 Modeled Equation');


