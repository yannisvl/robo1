%% *** Robot (kinematic) model parameters *** 
close all; 
l(1) = 10.0;  %% in cm 
l(2) = 10.0;  
l(3) = 10.0;

%% *** sampling period *** 
%% *** for the robot motion, kinematic simulation: 
dt = 0.001; %dt = 0.001; i.e. 1 msec)   

%% *** Create (or load from file) reference signals *** 
%% *** DESIRED MOTION PROFILE - TASK SPACE *** 
Tf=10.0; 	% 10sec duration of motion 
t=0:dt:Tf;  

%xd0,td0,yd1: initial/final end-point position --> desired task-space trajectory  
xd0 = 10.0;	
xd1 =  10.0; 
yd0 = 0; 
yd1 = 10.0;
zd0 = -20.0;
zd1 = 0;

% Example of desired trajectory : linear segment (x0,y0)-->(x1,y1); Time duration: Tf; 
disp('Initialising Desired Task-Space Trajectory (Motion Profile) ...'); %% 
disp(' ');   
xd(1) = xd0; 
yd(1) = yd0; 
zd(1) = zd0;

n = length(t);
vx = zeros(n,1);
vy = zeros(n,1);
vz = zeros(n,1);

vy(1) = 0; 
vz(1) = 0;
lambda_x = (xd1-xd0)/Tf; 
lambda_y = (yd1-yd0)/Tf;
lambda_z = (zd1 - zd0) / Tf;
kmax=Tf/dt + 1; 
for k=2:kmax   
   xd(k) = xd(k-1) + lambda_x*dt; 
   vx(k) = ( xd(k)-xd(k-1) ) / dt;
   yd(k) = yd(k-1) + lambda_y*dt;
   vy(k) = ( yd(k)-yd(k-1) ) / dt;
   zd(k) = zd(k-1) + lambda_z*dt;
   vz(k) = ( zd(k)-zd(k-1) ) / dt;
end  


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%% ****** KINEMATIC SIMULATION - Main loop ****** 
disp('Kinematic Simulation ...'); %% 
disp(' '); %%  

%% ***** INVESRE KINEMATICS  -->  DESIRED MOTION - JOINT SPACE ***** 
%% compute the reference joint-motion vectors: 
%% {qd(k,i), i=1,...,n (num of degrees of freedom), with k=1,..., kmax,} 
%% and reference joint (angular) velocities {qd_1(k,i)} 
rd2 = xd(:).^2 + yd(:).^2 + zd(:).^2; 
qd(:,3) = acos( (rd2(:)-l(2)^2-l(3)^2-l(1)^2)./(2*l(2)*l(3)) ); %% 1st solution: elbow down 	
%% or qd(:,2) = -acos( (rd2(:)-l(1)^2-l(2)^2)./(2*l(1)*l(2)) ); %% 2nd solution: elbow up 
s3 = sin(qd(:,3));
c3 = cos(qd(:,3));
qd(:,2) = 2*atan((l(2)+l(3)*c3 + sqrt(l(2)^2+l(3)^2+2*l(2)*l(3)*c3 - yd(:).^2)) ./ (l(3)*s3 + yd(:)));

s2 = sin(qd(:,2));
c2 = cos(qd(:,2));
c23 = cos(qd(:,2)+qd(:,3));
qd(:,1) = 2*atan((l(3)*c23 + l(2)*c2 + sqrt((l(2)*c23+l(2)*c2).^2+l(1)^2-yd(:).^2)) ./ (yd(:) + l(1)) );

%%inverse kinematics
s1 = sin(qd(:,1));
c1 = cos(qd(:,2));
qteleia(:,1) = (c1.*vx()+ vz.*s1) ./ (l(3)*c23+l(2)*c2)  ;

s23 = sin(qd(:,2)+qd(:,3));
a1 = l(3)*c23.*s1.*c23;
a2 = l(2).*s1.*c2.*c23;
a3 = l(1).*c1.*c23;
temp1 = (a1+a2+a3).*vx;
temp2 = c23.*(l(1).*s1 - l(3).*c1.*c23 - l(2).*c1.*c2).*vz;
temp3 = (l(3)*s3.*(l(3).*c23+l(2).*c2));
temp4 =  (s23.*vy) ./ (l(2)*s3);
qteleia(:,2) =((temp1 + temp2) ./ temp3) + temp4;

qteleia(:,3) = ((l(3)*s1.*c23+l(2)*s1.*c2+l(1)*c1).*vx + (l(3)*s23+l(2)*s2).*vy + (l(1)*s1 - l(3)*c1.*c23 - l(2)*c1.*c2).*vz)./ (-l(3)*l(2)*s3);
%% ***** FORWARD KINEMATICS  JOINT MOTION -->  CARTESIAN POSITIONS ***** 
%%(xd1, yd1, zd1) : cartesian position of the 2nd link's local reference frame 
xd1 = l(1)*cos(qd(:,1));   
yd1 = zeros(n,1);
zd1 = l(1)*sin(qd(:,1)); 
%%(xd2, yd2, zd1) : cartesian position of the 3d link's local reference frame 
xd2 = l(1)*cos( qd(:,1) ) + l(2)*cos( qd(:,1) ).*sin( qd(:,2) );   
yd2 = l(2)*sin ( qd(:,2) );
zd2 = l(1)*sin( qd(:,1) ) - l(2)*cos( qd(:,1) ).*cos( qd(:,2) );
%%(xd3, yd3) : cartesian position of the 3nd link's local reference frame
xd3 = l(3)*sin( qd(:,1) ).*cos( qd(:,2)+qd(:,3) ) + l(2)*sin( qd(:,1) ).*cos( qd(:,2) ) + l(1)*cos (qd(:,1) );
yd3 = l(3)*sin( qd(:,2)+qd(:,3) ) + l(2)*sin( qd(:,2) );
zd3 = -l(3)*cos( qd(:,1) ).*cos( qd(:,2)+qd(:,3) ) -  l(2)*cos( qd(:,1) ).*cos( qd(:,2) ) + l(1)*cos (qd(:,1) );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  


%% *** SAVE and PLOT output data *** %%** use functions plot(...)  
%%save; 
%% --> save data to 'matlab.mat' file   

fig1 = figure;  
subplot(3,3,1); 
plot(t,xd); 
ylabel('xd (cm)'); 
xlabel('time t (sec)');  

subplot(3,3,2); 
plot(t,yd); 
ylabel('yd (cm)'); 
xlabel('time t (sec)');  

subplot(3,3,3); 
plot(t,zd); 
ylabel('zd (cm)'); 
xlabel('time t (sec)');  

subplot(3,3,4); 
plot(t,qd(:,1)); 
ylabel('qd1 (rad)'); 
xlabel('time t (sec)');  

subplot(3,3,5); 
plot(t,qd(:,2)); 
ylabel('qd2 (rad)'); 
xlabel('time t (sec)');    

subplot(3,3,6); 
plot(t,qd(:,3)); 
ylabel('qd3 (rad)'); 
xlabel('time t (sec)'); 

subplot(3,3,7); 
plot(t,qteleia(:,1)); 
ylabel('qw1 (rad)'); 
xlabel('time t (sec)'); 

subplot(3,3,8); 
plot(t,qteleia(:,2)); 
ylabel('qw2 (rad)'); 
xlabel('time t (sec)'); 

subplot(3,3,9); 
plot(t,qteleia(:,3)); 
ylabel('qw3 (rad)'); 
xlabel('time t (sec)'); 
%%*** stick diagram --> animate robot motion ... (**optional**) 
%% within a for (or while) loop, use periodic plot(...) functions to draw the geometry (current pos)  
%% of the robot, and thus animate its motion ...  

fig2 = figure; 
axis([-5 25 -15 15]) %%set xyz plot axes (caution: square axes, i.e. dx=dy=dy) 
axis on 
hold on 
xlabel('x (cm)'); 
ylabel('y (cm)'); 
zlabel('z (cm)');
plot3(xd,yd,zd,'rs'); 
dtk=1000; %% plot robot position every dtk samples, to animate its motion 
plot3(0,0,0,'o'); 
hold on
for tk=1:dtk:kmax   %%% 
   pause(0.1);	%% pause motion to view successive robot configurations    			
   hold on
   plot3(xd1(tk),yd1(tk),zd1(tk),'o');    
   hold on 
   plot3([xd1(tk),xd2(tk)],[yd1(tk),yd2(tk)],[zd1(tk),zd2(tk)]);
   plot3(xd2(tk),yd2(tk),zd2(tk),'o'); 
   plot3([xd2(tk),xd3(tk)],[yd2(tk),yd3(tk)],[zd2(tk),zd3(tk)]);	
   plot3(xd3(tk),yd3(tk),zd3(tk),'o'); 
   plot3([xd3(tk),xd(tk)],[yd3(tk),yd(tk)],[zd3(tk),zd(tk)]);
   plot3(xd(tk), yd(tk), zd(tk), 'g+');  
end