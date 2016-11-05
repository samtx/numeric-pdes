Case=1;
clc
clear all
close all
S=100;
l=50;
q=200;
D=8.8E7;
A=(q*l^3)/(2*D);
gamma=((S*l^2)/D);
Case=1;
Num=10; %Num of elements
p=2; %Polynomial of degree
N=(Num*p)+1; %Nodes
h=1/Num;
Mass=(h/30)*[4 2 -1; 2 16 2; -1 2 4];
Stiff=(1/(3*h))*[7 -8 1; -8 16 -8; 1 -8 7];
t=linspace(0,1,N);
[Global_stiff, Global_mass] = Assembling_Matrices_Quadratic(Stiff,Mass,N,p,Case);  
t=linspace(0,1,N);
f=A.*t.*(1-t);
F=Global_mass*f';
A_matrix=Global_stiff+gamma*Global_mass
A_matrix(1,:)=[];
F(1)=[];
F(end)=[];
A_matrix(end,:)=[];
size(A_matrix)
size(f);
sol=A_matrix\F

uu=-t.^2+t-2/gamma;
vv=2/(gamma*sinh(sqrt(gamma)));
pp=sinh(sqrt(gamma).*t)+sinh(sqrt(gamma).*(1-t));
solution=A/gamma*(uu+vv*pp);
size(solution);
compare=[sol solution']
error=abs(compare(:,1)-compare(:,2));
H_one=error'*Global_stiff*error+error'*Global_mass*error;
L2=error'*Global_mass*error;
H_one
L2
plot(t,solution)
hold on
plot(t,sol,'r--')

