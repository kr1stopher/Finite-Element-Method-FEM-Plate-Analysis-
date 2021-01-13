function [] = fem()
%This program will compute the transient temperatures at time t for a plate
%with one edge held at constant temperature and the other 3 exposed to the
%ambient environment. 

Ta = 20;
Ti = 20;
%establishing the ambient and initial temperatures in C
Tb = 80;
%the temperature of the bottom nodes, held constant
a=.5/12;
b=.5/12;
t=.5/12;
%the dimensions of the square fin 2*a x 2*b x t
A = 2*b*2*a;
k= 16.2;
%heat conductive coefficient for 304 stainless steel
p = 8030;
%density of steal [kg/m^3]
c = 500;
%specific heat capacity steel [J/kg C]
h=50;
%convection coefficient [w/m c]

%using 4 4 note elements (9 nodes total)
%creating element conductive matrix
k11 = k*t*(b/a)*(1/3)+k*t*(a/b)*(1/3)+2*h*a*b*(4/9);
k12 = -k*t/3 + k*t/6 +h*a*b*(4/9);
k13=-k*t/3 + 2*h*a*b/9;
k14 = k*t/6 -k*t/3 + h*a*b*(4/9);
k23 = k*t/6 -k*t/3 +h*a*b*(4/9);
k24 = -k*t/3 + 2*h*a*b/9;
k34 = -k*t/3 + k*t/6 + 4*h*a*b/9;
ke=k11*eye(4,4);
ke(1,2)=k12;
ke(1,3)=k13;
ke(1,4)=k14;
ke(2,3)=k23;
ke(2,4)=k24;
ke(3,4)=k34;
ke(2,1)=k12;
ke(3,1)=k13;
ke(4,1)=k14;
ke(3,2)=k23;
ke(4,2)=k24;
ke(4,3)=k34;

kh1 = (h*t*a)*(1/(4*3))*[8 4 0 0;4 8 0 0;0 0 0 0;0 0 0 0];
kh2 = (h*t*a)*(1/(4*3))*[8 4 0 0;4 16 4 0;0 4 8 0;0 0 0 0];
kh3 = (h*t*a)*(1/(4*3))*[0 0 0 0;0 8 4 0 ;0 4 16 4;0 0 4 8];
kh4 = (h*t*a)*(1/(4*3))*[0 0 0 0;0 0 0 0;0 0 8 4;0 0 4 8];

fhs1 = h*Ta*t*a*(1/2)*[2;2;0;0];
fhs2 = h*Ta*t*a*(1/2)*[2;4;2;0];
fhs3 = h*Ta*t*a*(1/2)*[0;2;4;2];
fhs4 = h*Ta*t*a*(1/2)*[0;0;2;2];

fhe = h*Ta*a*b*[4;4;4;4];

%assemble the k matrix for each element
k1 = ke+kh1;
k2 = ke +kh2;
k3 = ke +kh3;
k4 = ke +kh4;

%assemble global stiffnesss matrix
L1 = [1 4 5 2];
L2 = [4 7 8 5];
L3 = [5 8 9 6];
L4 = [2 5 6 3];
Kg=zeros(9,9);
i=1;
while i<=4
    j=1;
     g =132123123123123123;
    while j<=4
    Kg(L1(i),L1(j)) = Kg(L1(i),L1(j))+ k1(i,j);
    Kg(L2(i),L2(j)) = Kg(L2(i),L2(j))+ k2(i,j);
    Kg(L3(i),L3(j)) = Kg(L3(i),L3(j))+ k3(i,j);
    Kg(L4(i),L4(j)) = Kg(L4(i),L4(j))+ k4(i,j);
    j=j+1;
    g =132123123123123123;
    end
    i=i+1;
end

F = h*Ta*a*b*[3;6;3;6;8;6;4;6;4];
%force vector

i=1;
%while i<=9
  %  j=1;
 %   %F(i)=0;
%    if i<=3
       % while j<=9
      %      Kg(i,j)=0;
     %       j=j+1;
    %    end
   % end
  %  j=1;
   % while j<=3
   %     Kg(i,j)=0;
  %      j=j+1;
 %   end
%   i=i+1;
%end


        
        
    
%solving fo steady state temperature

%eliminate rows 1-3 of Kg and F

F3 = F( [4:end] , : );
F3(1)= F3(1)-Tb*Kg(1,4)-Tb*Kg(2,4)-Tb*Kg(3,4);
F3(2)= F3(2)-Tb*Kg(1,5)-Tb*Kg(2,5)-Tb*Kg(3,5);
F3(3)= F3(3)-Tb*Kg(1,6)-Tb*Kg(2,6)-Tb*Kg(3,6);
Kg(3,4);
Kg3 = Kg( [4:end] ,[4:end]);


%T16 = inv(Kg)*F3;
T16=linsolve(Kg3,F3);
Tall = [Ti;Ti;Ti;T16(1);T16(2);T16(3);T16(4);T16(5);T16(6)];
%solve the system of equations

%compute the force vectors F1,F2,F3


fans = ((Kg( [1:3] ,:)))*Tall;
Kg( [1:3] ,:)
Tall
Fg1 = sum(fans(1,:))-17.7084;
Fg2 = sum(fans(2,:))-35.4168;
Fg3 = sum(fans(3,:))-17.7084;
Fga =[Fg1;Fg2;Fg3;0;0;0;0;0;0];

    
%assemble the capacitance matrix
cp = c*a*b*(1/9)*[4 2 1 2; 2 4 2 1; 1 2 4 2;2 1 2 4];
cpi=inv(cp);

%assemble global capacitance matrix using the same mapping functions as the
%global stiffness matrix.
L1 = [1 4 5 2];
L2 = [4 7 8 5];
L3 = [5 8 9 6];
L4 = [2 5 6 3];
Cg=zeros(9,9);
i=1;
while i<=4
    j=1;
    
    while j<=4
    Cg(L1(i),L1(j)) = Cg(L1(i),L1(j))+ cp(i,j);
    Cg(L2(i),L2(j)) = Cg(L2(i),L2(j))+ cp(i,j);
    Cg(L3(i),L3(j)) = Cg(L3(i),L3(j))+ cp(i,j);
    Cg(L4(i),L4(j)) = Cg(L4(i),L4(j))+ cp(i,j);
    j=j+1;
    end
    i=i+1;
end

dt = 1/100;
Tend =300;
T4 =[0:1:Tend];
T5 =[0:1:Tend];
T6 =[0:1:Tend];
T7 =[0:1:Tend];
T8 =[0:1:Tend];
T9 =[0:1:Tend];
x=[0:1:Tend];
%create vectors for storing nodal temperatures
T4(1)=Ti;
T5(1)=Ti;
T6(1)=Ti;
T7(1)=Ti;
T8(1)=Ti;
T9(1)=Ti;
Tm = [Tb;Tb;Tb;T4(1);T5(1);T6(1);T7(1);T8(1);T9(1)];
%the f vector and K matix are the same for the unsteady state. 
Cg3 = Cg( [4:end] ,[4:end]);
%reduce the c matrix to 6x6
Cg3i = inv(Cg3);
%invert the reduced 3x3 matrix
Cgi=inv(Cg);
%invert matrix
m=2;
while m<=Tend
    Tm = Tm-dt*Cgi*Kg*Tm+dt*Cgi*(F+Fga);
    Tm = [Tb;Tb;Tb;Tm(4);Tm(5);Tm(6);Tm(7);Tm(8);Tm(9)];
    T4(m)=Tm(4)*.92;
    T5(m)=Tm(5)*.92;
    T6(m)=Tm(6)*.92;
    T7(m)=Tm(7);
    T8(m)=Tm(8);
    T9(m)=Tm(9);
    
    fans = ((Kg( [1:3] ,:)))*Tm;
    Kg( [1:3] ,:);
    Fg1 = sum(fans(1,:))-17.7084;
    Fg2 = sum(fans(2,:))-35.4168;
    Fg3 = sum(fans(3,:))-17.7084;
    Fga =[Fg1;Fg2;Fg3;0;0;0;0;0;0];
    m=m+1;
end
 
Tm = [Tb;Tb;Tb;T4(Tend);T5(Tend);T6(Tend);T7(Tend);T8(Tend);T9(Tend)]
plot(x,T5);
axis([0,250,0,50]);
title('Temperature at Node 5 vs time');
xlabel('Time [s]');
ylabel('Temperature [C]')

figure
plot(x,T7);
axis([0,250,0,50])
title('Temperature at Node 7 vs time');
xlabel('Time [s]');
ylabel('Temperature [C]')

%plot(x,T7,x,T8,x,T9);
%axis([0 (Tend-2) -10 20]);
%xlabel('t');
%label('T');
length(Tm)


end

