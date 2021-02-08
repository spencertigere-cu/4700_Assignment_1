%Spencer Tigere
%101001717

clc
clear
clearvars
set(0,'DefaultFigureWindowStyle','docked')

%Constants
temp = 300;
c = 299792458;
force_grav = 9.80665;
H = 1.054571596e-34;
H = H * 2 * pi;
Q  = 1.60217653e-19;
m_0 = 9.10938215e-31;
m_eff = 0.26 * m_0;
boltz_c = 1.3806504e-23;
permittivity = 8.854187817e-12;
mu = 1.2566370614e-6;
am = 1.66053892e-27;

%Particle movement. Question 1
thermal_vel = sqrt((boltz_c * 300) / m_eff); 
stand_dev = thermal_vel/(sqrt(2));
delta_t = 7.5 * 10 ^ -15;
array_1 = zeros(1, 1000);
tmp_arary_1 = (1:1:1000);
s_m = 300;
w_x = 200 * 10 ^ -9;
l_x = 100 * 10 ^ -9;
size_particles = 50;
tmn = 0.2 * 10 ^ -12;

%current x and y positions
x_curr = rand(1, size_particles) .* w_x;
y_curr = rand(1, size_particles) .* l_x;

theta = rand(1,size_particles)*2*pi;

%x and y velocities
x_vel = randn(1,size_particles) .*stand_dev;
y_vel = randn(1,size_particles) .*stand_dev;
 

for i = 1:1000
    
    x_1 = x_curr;
    x_curr = x_1 + delta_t*x_vel;
    y_1 = y_curr;
    y_curr = y_1 + delta_t*y_vel;
    j1 = (y_curr >= l_x);
    j2 = (y_curr < 0);
    
    y_vel(j1) = -y_vel(j1);
    y_vel(j2) = -y_vel(j2);
    j3 = ((x_curr) >= w_x);
    j4 = ((x_curr) < 0); 
    
    x_1(j3) = x_1(j3) - w_x;
    x_curr(j3) = x_curr(j3) - w_x;
    x_1(j4) = x_1(j4) + w_x;
    x_curr(j4) = x_curr(j4) + w_x;
    inbox = ((x_curr <= (1.30 * w_x/2) & (x_curr >= (0.90 * w_x/2))) & ((y_curr < (l_x/3)) | y_curr >= (2*l_x/3)));
    middle = ((x_1<= (1.30 * w_x/2) & (x_1>= (0.90 * w_x/2))) & ((y_1 > (l_x/3)) & y_1 <= (2*l_x/3)));
    
    plot(x_curr,y_curr,'r.');
    pause(0.05)
    hold on
    
    line([0.90*w_x/2 0.90*w_x/2], [l_x 2*l_x/3]);
    line([1.30*w_x/2 1.30*w_x/2], [l_x 2*l_x/3]);
    line([0.90*w_x/2 1.30*w_x/2], [l_x l_x]);
    line([0.90*w_x/2 1.30*w_x/2], [2*l_x/3 2*l_x/3]);
    
    line([0.90*w_x/2 0.90*w_x/2], [0 l_x/3]);
    line([1.30*w_x/2 1.30*w_x/2], [0 l_x/3]);
    line([0.90*w_x/2 1.30*w_x/2], [0 0]);
    line([0.90*w_x/2 1.30*w_x/2], [l_x/3 l_x/3]);
    
    vrms = sqrt ((x_vel.^2)+(y_vel.^2));
    x_vel(inbox&(~middle)) = -x_vel(inbox&(~middle));
    y_vel(inbox&middle) = -y_vel(inbox&middle);
    meanfreetime = (thermal_vel * Q)/200;
    
    %Question 2
    meanfreepath = mean(vrms) * meanfreetime;
    vrms = sqrt ((x_vel.^2)+(y_vel.^2));
    s_m = (sqrt (2) * (mean(vrms)^2) * m_eff )/boltz_c;
    array_1 (1,i)=  s_m;
    
    [x_mesh, y_mesh] = meshgrid(0:(w_x/10):w_x, 0:(l_x/10):l_x);
    electron_mat = zeros(11, 11);
    temperture_mat = zeros(11, 11);
    numelec_t = 0;
    t_vel = 0;
    
    for j = 1:10
        efx_min = x_mesh(1, j);
        efx_max = x_mesh(1,j+1);
        for k = 1:10
            efy_min = y_mesh(k, 1);
            efy_max = y_mesh(k+1, 1);
            for m = 1:size_particles
                if((x_curr(m) > efx_min) && (x_curr(m) < efx_max) && ((y_curr(m) > efy_min) && y_curr(m) < efy_max))
                    numelec_t = numelec_t + 1;
                    electron_mat(j, k) = electron_mat(j, j) + 1;
                    t_vel = t_vel + sqrt((x_vel(m) .^ 2) + (x_vel(m) .^ 2));
                    temperture_mat(j, k) = ((sqrt(2)*(t_vel/numelec_t) ^ 2) * m_eff) / boltz_c;
                end
            end
            t_vel = 0;
            numelec_t = 0;
        end
    end
end
 
fprintf("Mean Free Time = %f\n", meanfreetime);
fprintf("Mean Free Path = %f\n", meanfreepath);

figure (1)
title(["Average Temperature Value = " num2str(s_m)]);
xlabel("X-axis");
ylabel("Y-axis");

figure (2)
plot(tmp_arary_1, array_1);
title('Temperature Across Time');
ylabel('Temperature');
xlabel('Time');
hold on

figure(3); Histogram(vrms, 15);
title('Histogram of Thermal Velocities');
xlabel("X-Axis");
ylabel("Y-Axis");

figure(4); surf(electron_mat);
title('Density Mapping');
xlabel("X-Axis");
ylabel("Y-Axis");

figure(5); surf(temperture_mat);
title('Temperature Mapping');
xlabel("X-Axis");
ylabel("Y-Axis");
 
clc
clear
clearvars
set(0,'DefaultFigureWindowStyle','docked')
particles = 1000;
 
global C X Y
C.Q_0 = 1.60217653e-19;             
C.Hb = 1.054571596e-34;             
C.H = C.Hb * 2 * pi;                
C.m_0 = 9.10938215e-31;             
C.kb = 1.3806504e-23;               
C.eps_0 = 8.854187817e-12;          
C.mu_0 = 1.2566370614e-6;           
C.c = 299792458;                    
C.g = 9.80665;                      
C.m_n = 0.26*C.m_0;                 
 
region_x = 200e-9;
region_y = 100e-9;
step_size = 1e-9;
timestep = 1000;
T = 300;
v_tH = sqrt(2*C.kb*T/C.m_n);
v_cHange = step_size/v_tH;
MT_C = 0.2e-12;
 
X = rand(2,particles);
Y = rand(1,particles);
x_currition(1,:) = X(1,:)*region_x;
y_currition(1,:) = Y(1,:)*region_y;
 
cHeck_X_left = x_currition > 0.8e-7;
cHeck_X_rigHt = x_currition < 1.2e-7;
cHeck_X = cHeck_X_left & cHeck_X_rigHt;
cHeck_top = y_currition > 0.6e-7;
cHeck_bottom = y_currition < 0.4e-7;
box_top = cHeck_top & cHeck_X;
box_bottom = cHeck_bottom & cHeck_X;
IN_A_BOX = box_top | box_bottom;
 
while(sum(IN_A_BOX) > 0)
    
    temp_x = rand(1,sum(IN_A_BOX));
    temp_y = rand(1,sum(IN_A_BOX));
    x_currition(IN_A_BOX) = temp_x*region_x;
    y_currition(IN_A_BOX) = temp_y*region_y;
    
    cHeck_X_left = x_currition > 0.8e-7;
    cHeck_X_rigHt = x_currition < 1.2e-7;
    cHeck_X = cHeck_X_left & cHeck_X_rigHt;
    cHeck_top = y_currition > 0.6e-7;
    cHeck_bottom = y_currition < 0.4e-7;
    box_top = cHeck_top & cHeck_X;
    box_bottom = cHeck_bottom & cHeck_X;
    IN_A_BOX = box_top | box_bottom;
end

angle(1,:) = X(2,:)*2*pi;
sigma = sqrt(C.kb*T/C.m_n)/4;
max_boltz_dist = makedist('Normal',v_tH,sigma);
velocity = random(max_boltz_dist,1,particles);

figure(6)
Histogram(velocity)
title('Particle Velocity Histogram')
X_velocity = v_cHange*velocity(1,:).*cos(angle(1,:));
Y_velocity = v_cHange*velocity(1,:).*sin(angle(1,:));
PSCAT = 1 - exp(-v_cHange/MT_C);
mfp_vec = zeros(1,particles);
clc
clear
clearvars
set(0,'DefaultFigureWindowStyle','docked')

K = 1.3806e-23;
m = 0.26*9.1093e-31;
 
x = 200e-9*rand(1000,1);
y = 100e-9*rand(1000,1);
 
yboundSpecular = true; 
xboundSpecular = true;
boxSpecular = true;
 
inbox1 = x > 80e-9 & x < 120e-9 & y > 60e-9;
inbox2 = x > 80e-9 & x < 120e-9 & y < 40e-9;
x(inbox1) = x(inbox1) + ((rand() > 0.5)*2 - 1)*40e-9*rand(size(x(inbox1)));
x(inbox2) = x(inbox2) + ((rand() > 0.5)*2 - 1)*40e-9*rand(size(x(inbox2)));
y(inbox1) = y(inbox1) - 0.2*rand(size(y(inbox1)));
y(inbox2) = y(inbox2) + 0.2*rand(size(y(inbox2)));
 
T = 300;
vtH = sqrt(2*K*T/m);
std = sqrt(K*T/m);
Vx = normrnd(0,std,[1000,1]);
Vy = normrnd(0,std,[1000,1]);
V = sqrt(Vx.^2 + Vy.^2);
Tplot = zeros(1000,1);
figure(1)
Histogram(V);
dt = 0.5e-14;
 
xold = x;
yold = y;
for i =1:400
    %Defines tHe boundaries of tHe simulation as well as tHe boxes
    yboundTop = y > 100e-9;
    yboundBottom = y < 0;
    inbox1 = x >= 80e-9 & x <= 120e-9 & y >= 60e-9;
    inbox2 = x >= 80e-9 & x <= 120e-9 & y <= 40e-9;
    xboundRigHt = x > 200e-9;
    xboundLeft = x < 0;

    if xboundSpecular
        Vx(xboundRigHt | xboundLeft) = - Vx(xboundRigHt | xboundLeft);
    else
        theta = pi*rand();
        Vx(xboundRigHt | xboundLeft) = V(xboundRigHt | xboundLeft)*cos(theta);
        Vy(xboundRigHt | xboundLeft) = V(xboundRigHt | xboundLeft)*sin(theta);
    end
    
    %Reflection off of y boundary
    if yboundSpecular
        Vy(yboundTop | yboundBottom) = -Vy(yboundTop | yboundBottom);
    else
        theta = pi*rand();
        Vy(yboundTop | yboundBottom) = V(yboundTop | yboundBottom)*cos(theta);
        Vx(yboundTop | yboundBottom) = V(yboundTop | yboundBottom)*sin(theta);
    end

    if boxSpecular
        Vx(inbox1 & yold >= 60e-9) = -Vx(inbox1 & yold >= 60e-9);
        Vx(inbox2 & yold <= 40e-9) = -Vx(inbox2 & yold <= 40e-9);
        
        Vy(inbox1 & yold <= 60e-9) = -Vy(inbox1 & yold <= 60e-9);
        Vy(inbox2 & yold >= 40e-9) = -Vy(inbox2 & yold >= 40e-9);
    else
        theta = pi*rand();
        
        Vx(inbox1 & yold >= 60e-9) = V(inbox1 & yold >= 60e-9)*cos(theta);
        Vx(inbox2 & yold <= 40e-9) = V(inbox2 & yold <= 40e-9)*cos(theta);
        Vy(inbox1 & yold >= 60e-9) = V(inbox1 & yold >= 60e-9)*sin(theta);
        Vy(inbox2 & yold <= 40e-9) = V(inbox2 & yold <= 40e-9)*sin(theta);

        Vy(inbox1 & yold <= 60e-9) = V(inbox1 & yold <= 60e-9)*cos(theta);
        Vy(inbox2 & yold >= 40e-9) = V(inbox2 & yold >= 40e-9)*cos(theta);
        Vx(inbox1 & yold <= 60e-9) = V(inbox1 & yold <= 60e-9)*sin(theta);
        Vx(inbox2 & yold >= 40e-9) = V(inbox2 & yold >= 40e-9)*sin(theta);
    end
    y(yboundTop) = 100e-9;
    y(yboundBottom) = 0;
    x(xboundRigHt) = 200e-9;
    x(xboundLeft) = 0;
    x(inbox1 & yold >= 60e-9 & x <= 100e-9) = 80e-9;
    x(inbox1 & yold >= 60e-9 & x > 100e-9) = 120e-9;
    x(inbox2 & yold <= 40e-9 & x <= 100e-9) =80e-9;
    x(inbox2 & yold <= 40e-9 & x >= 100e-9) =120e-9;
    y(inbox1 & yold <= 60e-9) = 60e-9;
    y(inbox2 & yold >= 60e-9) = 40e-9;
   
 
    xold = x;
    yold = y;
    x = x + Vx*dt;
    y = y + Vy*dt;
    scatter = rand(1000,1) < (1 - exp(-dt/0.2e-12));
    Vx(scatter) = normrnd(0,std,size(Vx(scatter)));
    Vy(scatter) = normrnd(0,std,size(Vy(scatter)));
    xplot = transpose([xold(1:20) x(1:20)]);
    yplot = transpose([yold(1:20) y(1:20)]);
    Tplot(i) = (1/(2*K))*mean(Vx.^2 + Vy.^2)*m;
    xlim([0 200e-9])
    ylim([0 100e-9])
    drawnow  
end

temp_sum_x = zeros(20,10);
temp_sum_y = zeros(20,10);
temp_num = zeros(20,10);
 
for i=1:1000
 x1 = floor(x(i)/1e-8);
 y1 = floor(y(i)/1e-8);
 if(x1<=0)
 x1 = 1;
 end
 if(y1<=0)
 y1= 1;
 end
 if(y1>100)
     y1 = 100;
 end
 if(x1>200)
     x1=200;
 end
 temp_sum_y(x1,y1) = temp_sum_y(x1,y1) + Vy(i).^2;
 temp_sum_x(x1,y1) = temp_sum_x(x1,y1) + Vx(i).^2;
 temp_num(x1,y1) = temp_num(x1,y1) + 1;
 
end
 
temp = (temp_sum_x + temp_sum_y).*m./K./2./temp_num;
temp(isnan(temp)) = 0;
temp = transpose(temp);
 
 
figure(7)
surf(temp)
title('Temperature Map');
xlabel('X-axis(nm)');
ylabel('Y-axis(nm)');