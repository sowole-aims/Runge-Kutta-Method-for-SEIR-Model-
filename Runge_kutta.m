% It calculates SEIR Model using Runge-Kutta 4th order method
% Authour by: Sowole Oladimeji Samuel, African Institute for Mathematical Sciences, Senegal.
% Date 11/06/2018.
clc;                                   % Clears the screen
clear all;

%%%% specifying Parameters used %%%%
 b = 0.0534; beta = 0.09091; sigma = 0.00; mu =  0.03835; gamma = 0.125;  alpha = 0.14286;                              % step size
h= 0.0001;  tfinal = 10;
N = ceil(tfinal/h)                           % Calculates upto y(10000) ceil(tfinal/h)

                               % initial condition
f1 = @(t,S, E, I,R) b - (beta*I + (mu))* S;      % change the function as you desire
f2 = @ (t, S,E,I,R)  (beta * S * I ) - (mu + alpha + sigma  )*E
f3 = @ (t, S,E,I,R)  (alpha * E)  - (mu + gamma )*I
f4 = @ (t, S,E,I,R)  (gamma * I) + (sigma*E) - (mu*R)  
% initial conditions
t(1)=0; S(1) = S0 = 800000;  E(1)= E0 = 150000; I(1)=I0=10000; R(1) =R0 =40000;
%step size
for i=1:N;  % update loop
%update time
t(i+1) = t(i) + h; 

 %%%%% Stage One %%%%%
    
 K1S = f1( t(i) ,  S(i),  E(i),  I(i),  R(i)   );
 K1E = f2( t(i) ,  S(i) ,  E(i) , I(i) ,  R(i) );
 K1I = f3( t(i) , S(i) , E(i) , I(i) , R(i) );
 K1R = f4( t(i) , S(i) , E(i) , I(i) , R(i) );
 %%%%% Stage Two %%%%%%
 K2S = f1( t(i) + 0.5*h,   S(i) + 0.5*h*K1S,   E(i) +  0.5*h*K1E,   I(i) + 0.5*h* K1I ,   R(i) + 0.5*h* K1R );
 K2E = f2( t(i) + 0.5*h,   S(i) + 0.5*h*K1S,   E(i) +  0.5*h*K1E,   I(i) + 0.5*h* K1I ,   R(i) + 0.5*h* K1R );
 K2I = f3( t(i) + 0.5*h,   S(i) + 0.5*h*K1S,   E(i) +  0.5*h*K1E,   I(i) + 0.5*h* K1I ,   R(i) + 0.5*h* K1R );
 K2R = f4( t(i) + 0.5*h,   S(i) + 0.5*h*K1S,   E(i) +  0.5*h*K1E,   I(i) + 0.5*h* K1I ,   R(i) + 0.5*h* K1R );

     %%%%% Stage Three %%%%%
 K3S = f1( t(i) + 0.5*h,   S(i) + 0.5*h*K2S,   E(i) +  0.5*h*K2E,   I(i) + 0.5*h* K2I ,   R(i) + 0.5*h* K2R );
 K3E = f2( t(i) + 0.5*h,   S(i) + 0.5*h*K2S,   E(i) +  0.5*h*K2E,   I(i) + 0.5*h* K2I ,   R(i) + 0.5*h* K2R );
 K3I = f3( t(i) + 0.5*h,   S(i) + 0.5*h*K2S,   E(i) +  0.5*h*K2E,   I(i) + 0.5*h* K2I ,   R(i) + 0.5*h* K2R );
 K3R = f4( t(i) + 0.5*h,   S(i) + 0.5*h*K2S,   E(i) +  0.5*h*K2E,   I(i) + 0.5*h* K2I ,   R(i) + 0.5*h* K2R );
 %%%%%  Stage Four %%%%%
 K4S = f1( t(i) + h,   S(i) + h*K3S,   E(i) +  h*K3E,   I(i) + h* K3I ,   R(i) + h*K3R );
 K4E = f2( t(i) + h,   S(i) + h*K3S,   E(i) +  h*K3E,   I(i) + h* K3I ,   R(i) +  h*K3R );  
 K4I = f3( t(i) + h,   S(i) + h*K3S,   E(i) +  h*K3E,   I(i) + h* K3I ,   R(i) +  h*K3R ); 
 K4R = f4( t(i) + h,   S(i) + h*K3S,   E(i) +  h*K3E,   I(i) +h* K3I ,   R(i) +  h*K3R );   
    
     %%%%% Now, the main equations %%%%%
     
    S(i+1) = S(i) + (1/6)*(K1S+ 2*K2S+2*K3S+K4S )*h;  
    E(i+1)  = E(i) + (1/6)*(K1E+2*K2E+2*K3E+K4E)*h; 
    I(i+1) = I(i) + (1/6)*(K1I+2*K2I+2*K3I+K4I)*h;
    R(i+1) = R(i) + (1/6)*(K1R+2*K2R+2*K3R+K4R)*h;
    %S(i+1)
    %E(i+1)
 
    
end

%plot the solutions
%figure(1); clf(1)
%plot(t,S)

%hold on
%plot(t,E)

%plot(t,I)
%plot(t,R)


%legend('S(t)')
%legend('E(t)')
%legend('I(t)')
%legend('R(t)')

%  The function to plots all together
plot(t,S,t,E,t,I,t,R)
legend('S(t)', 'E(t)', 'I(t)', 'R(t)')
xlabel('Time(years)')
ylabel('Populations')
set(gca, 'Fontsize', 12)
