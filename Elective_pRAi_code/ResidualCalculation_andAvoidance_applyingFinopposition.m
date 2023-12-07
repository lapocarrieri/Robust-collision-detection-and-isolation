clc
clear all
close all
load('InitializationSpatial_3R_robot')
G(3)=subs(G(3),{a(8)},{9.810000000000000});
G(2)=subs(G(2),{a(7),a(8)},{49.050000000000004,9.810000000000000});
stop=false;
SS=vpa(subs(C,{a(1),a(2),a(3),a(4),a(5),a(6),a(7),a(8)},{0.375000000000000,2.058333333333334,0.254166666666667,0.500000000000000,2.387500000000000,0.279166666666667,49.050000000000004,9.810000000000000}));
MM=vpa(subs(M,{a(1),a(2),a(3),a(4),a(5),a(6),a(7),a(8)},{0.375000000000000,2.058333333333334,0.254166666666667,0.500000000000000,2.387500000000000,0.279166666666667,49.050000000000004,9.810000000000000}));
JJ=vpa(subs(JJ,{L1,L2,L3},{0.5,0.5,0.4}));
Q=[1 0 1]';
dq0=[pi/8,pi/7,-pi/6];
sigma = zeros(3,1);
dQ=dq0';
index=1;
t0=0;
r=[0 0 0]'
tf=5;
DeltaT=0.001;
friction=[0 0 0]';
acc=[0 0 0]';
accs=[];
samples1=50;
samples2=50;
TauExternalForce=[1 2 5]';
n=3000;
cells=cell(1,n);
matrix = zeros(n, 3);
 S1=40; T1=sqrt(S1);
 S2=116 ; T2=sqrt(S2)*2;

gain=40*diag([1,1,1]);
gainE=100;
GainInv=inv((eye(3)+gain*DeltaT))*gain;
pinvJ=vpa(pinv(JJ(1:3,:)));
Q_sampled=matrix;
QD_sampled=matrix;
QDD_sampled=matrix;
TAU_applied=matrix;
ExternalTauCalculated=matrix;
ExternalTauCalculated2=matrix;
TauExternal=matrix;
time=zeros(n,1);

B_sampled=cells;
S_sampled=cells;
g_sampled=cells;
BQD=cells;
h_sampled=cells;
omega=pi;
q0=[1 0 1]';
F=[0 0 0]';
while (t0<tf)%(frequency * (t0) < 2*pi) % it ends when a circle is completed
     disp('time instant:')
     disp(t0);
     for i=1:3
        dQ(i)=dQ(i)+acc(i)*DeltaT;
     end
     for i=1:3
        Q(i)=Q(i)+dQ(i)*DeltaT;
     end
%     if (t0&1)==0
%         gain=gain*100;
%         S1=S1*4; T1=sqrt(S1);
%     end

%      TauExternalForce=[0 0 0]';
%      if t0>0.5
%         TauExternalForce=[1 4 8]';
%      end
     J=vpa(subs(JJ(1:3,:),{q1,q2,q3},{Q(1),Q(2),Q(3)}));
     if t0>0.2
        F=[1.5 -2 2.5]';
     end
    TauExternalForce=J'*F;
        g=subs(G,{q1,q2,q3},{Q(1),Q(2),Q(3)});
     M=double(subs(MM,{q1,q2,q3},{Q(1),Q(2),Q(3)}));
     S=double(subs(SS,{q1,q2,q3,dq1,dq2,dq3},{Q(1),Q(2),Q(3),dQ(1),dQ(2),dQ(3)}));
     q0=[cos(omega*t0);sin(omega*t0);1];
     pinvJ=vpa(subs(pinvJ,{q1,q2,q3},{Q(1),Q(2),Q(3)}));
     err_dp=dQ-dq0';
     err_p=Q-q0;
     Kd=0.1;
     Kp=1;
     Kr=0; % >0 if we want to apply a counterreaction to the force applied

     d2p_ref =  Kd * err_dp+Kp*err_p;  
     Uref_task = double((pinvJ * (d2p_ref)));
     
     TauFL = double(g + S*dQ + M * Uref_task);  % this is the applied torque 
%      ExternalForceAppliedActualFrame=[3 4 5]';
%      point=[0.1 0.2 -0.3];
%     S_fext =[0 -ExternalForceAppliedActualFrame(3) ExternalForceAppliedActualFrame(2) ; ExternalForceAppliedActualFrame(3) 0 -ExternalForceAppliedActualFrame(1) ; -ExternalForceAppliedActualFrame(2) ExternalForceAppliedActualFrame(1) 0 ];
%     m=-S_fext*point(1:3);
%     J_withwrenches=
%    TauExternalForce=(J_withwrenches'*[ExternalForceAppliedActualFrame;m])';
    
    Tau = TauFL+TauExternalForce -Kr*r;
    acc = M\ double(Tau - friction - S*dQ - g)  ;
      
    



    %PointOfApplication(index,:)=QtoP(Q_sampled(end,:),link)*point;
    
    noise =0.00001*[1 1 1]'* randn();
    %noise=0;
    Q_sampled(index,:)=Q;
    QD_sampled(index,:)=dQ+noise;
    QDD_sampled(index,:)=double(acc);
    TAU_applied(index,:)=(TauFL');
    B_sampled{index}=double(M);
    S_sampled{index}=double(S);
    g_sampled{index}=double(g);
    h_sampled{index}=double(S'*dQ-g);
    BQD{index}=double(M*QD_sampled(index, :)');
        if index==0
         fprintf('Starting of the Residual calculation:\n')

        end
     
         if index>samples1 % in order to have enough samples to calculate the residuals
            
            
             
            dsigma=[0 0 0]';
            sumTau=0;
            sumH=0;
            sumRes=0;
            p0=0;

            sumEdot=0;
            
             sumSigma=0;
         
%% Sliding Mode Observer
            for tt = 1:samples2+1
           
            t=index-samples2+tt-1;
            p=B_sampled{t}*QD_sampled(t, :)';
                if tt == 1
                   
                   p_hat = p;
                   %sigma = zeros(3,1);
                   
                else
                    
                    p=p_hat-p;
                    
                   signP=tanh(p*50);
%                  
%                    dp_hat=vpa(TAU_applied(t,:)'+h_sampled{t}+sigma-T2*p-T1*signP,3);
%                    
%                    dsigma=vpa(dsigma-S1*signP-S2*p,3);
                   dp_hat=double(TAU_applied(t,:)'+h_sampled{t}+sigma-T2*p-T1*signP);
                   
                   dsigma=double(-S1*signP-S2*p);
                   
                 
                   p_hat=p_hat+DeltaT*dp_hat;
                  
                   sigma=sigma+DeltaT*dsigma;
                   
                    
               end
            end
            
            ExternalTauCalculated(index,:)=vpa(sigma',3);
            error1=norm(sigma-TauExternalForce)
            ForceCalculated1(index,:)=pinv(J')*sigma;
            error2=norm(F-ForceCalculated1(index,:))
       
   

%% Momentum-based isolation of collisions
        

        for tt = 1:samples1+1

            t=index-samples1+tt-1;
          
            h=h_sampled{t};
            
            sumTau = sumTau + TAU_applied(t,:)';
            sumH = sumH + h;

        
                if tt == 1 
                   
                   p0 = B_sampled{t}*QD_sampled(t, :)';
                   
                   r = zeros(3,1);
                   
               else
                   %r = (eye(3)+gain*DeltaT) \( gain * ((B_sampled{t}*QD_sampled(t, :)' - p0) - (sumTau + sumH + sumRes)*DeltaT));
                   r = GainInv* ((BQD{t} - p0) - (sumTau + sumH + sumRes)*DeltaT);

                   sumRes = sumRes + r;
                   
                   
                   
               end
        end
        

        ExternalTauCalculated2(index,:)=r';
        ForceCalculated2(index,:)=pinv(J')*r;
        
        error3=norm(r-TauExternalForce)
        error4=norm(F-ForceCalculated2(index,:))
         end
        index=index+1;
        t0=t0+DeltaT;
        time(index)=t0;
        TauExternal(index,:)=TauExternalForce;
        
end
figure()
plot(time(samples2:index-1), TauExternal(samples2:index-1,:)', 'r', 'LineWidth', 2);
hold on;
plot(time(samples2+1:index-1), ExternalTauCalculated(samples2+1:index-1,:), 'g', 'LineWidth', 0.1);
hold off

figure()
plot(time(samples2+1:index-1), ForceCalculated1(samples2+1:index-1,:)', 'r', 'LineWidth', 2);
hold on;
plot(time(samples2+1:index-1), ForceCalculated2(samples2+1:index-1,:), 'b', 'LineWidth', 0.1);
hold off
plot(time(samples2+1:index-1), ForceCalculated2(samples2+1:index-1,:), 'b', 'LineWidth', 0.1);
hold off

figure();
scatter3(Q_sampled(:, 1), Q_sampled(:, 2), Q_sampled(:, 3), 'b.');
grid on;
xlabel('X-axis');
ylabel('Y-axis');
zlabel('Z-axis');
title('3D Scatter Plot');

hold off
figure()
plot(time(samples2:index-1), QD_sampled(samples2:index-1,:)', 'r', 'LineWidth', 2);


figure()
plot(time(samples2:index), TauExternal(samples2:index,:)', 'r', 'LineWidth', 2);
hold on;
plot(time(samples2+1:index-1), ExternalTauCalculated2(samples2+1:index-1,:), 'b', 'LineWidth', 0.1);