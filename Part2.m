%ELEC 4700 - Assignment 1 Resubmission
%Spencer Tigere 101001717
clc
clear all
clearvars

% Simulation variables    
electron_mass = 9.10938215e-31;    
effective_mass = 0.26*electron_mass;
dirac_constant = 1.054571596e-34;            
plank_constant = dirac_constant * 2 * pi;                   
boltzmann_constant = 1.3806504e-23;

r_w = 2e-7;% width of the region
r_l = 1e-7;% length of the region
T = 300;% Temperature (K)

thermal_vel = sqrt((2*boltzmann_constant*T)/effective_mass); % Thermal velocity
%I lost marks because I didnt have MFP or mean time. This has been
%resolved.
tmn = 2e-13;% Mean time between collisions
freepath = thermal_vel*2e-13;
fprintf('Using a temperature of 300 K, the thermal velocity is %f km/s. \n', thermal_vel/1000)
fprintf('Using a mean time to be 0.2 ps, the mean free path is %f nm. \n', freepath*1e9)
num_E = 10;
E_xpos = rand(1,num_E).*r_w;
E_ypos = rand(1,num_E).*r_l;
angle = rand(1,num_E).*2*pi;

%I lost marks here in my first attempt because I was 
%missing a histogram. that has been resolved
sig = sqrt(boltzmann_constant*T/effective_mass)/4;
MBdist = makedist('Normal',thermal_vel,sig);
E_vel = random(MBdist,1,num_E);
figure(1)
grid on
title('Particle Velocity Histogram')
hist(E_vel)
angle = rand(1,num_E).*2*pi;


E_xvel = E_vel.*cos(angle);
E_yvel = E_vel.*sin(angle);
d_t = 1e-9/thermal_vel;

p_scat = 1 - exp(-d_t/tmn);
E_vel = sqrt(E_xvel.^2 + E_yvel.^2);
counter = 0;

for i=1:1000
      for j = 1:num_E 
            if p_scat > rand
                counter = counter + 1;
                angleNew = rand(1).*2*pi;
                E_xvel(1,j) = random(MBdist,1).*cos(angleNew);
                E_yvel(1,j) = random(MBdist,1).*sin(angleNew);
       
            end
            
            if E_ypos(1,j) + E_yvel(1,j).*d_t  >= r_l || E_ypos(1,j) + E_yvel(1,j).*d_t <= 0
               E_yvel(1,j) = E_yvel(1,j)*-1;  
            end
            ypos(i,j) = E_ypos(1,j) + E_yvel(1,j).*d_t;
            if E_xpos(1,j) + E_xvel(1,j)*d_t  >= r_w || E_xpos(1,j) + E_xvel(1,j)*d_t <= 0
                if E_xpos(1,j) + E_xvel(1,j)*d_t >= r_w
                    xpos(i,j) = 0;
                else
                    xpos(i,j) = r_w;
                end
            end
                xpos(i,j) = E_xpos(1,j) + E_xvel(1,j).*d_t;
        end
    
    E_xpos = xpos(i,:);
    E_ypos = ypos(i,:);
    
    t1 = sqrt(E_xvel(1,:).^2 + E_yvel(1,:).^2);
    temp(i) = ((mean(t1)^2)*effective_mass)/(2*boltzmann_constant);% Semiconductor temperature
    
end
mfp = (1000/counter)*d_t*mean(E_vel);
meantime = d_t*(1000/counter);

%Electron Trajectory
%I lost marks here in my first submission because there were boxes in my
%simulation. That issue has been resolved here
figure (2)
plot(xpos,ypos,'-','LineWidth',2)
xlim([-0.1e-7 2.1e-7])
ylim([-0.1e-7 1.1e-7])
grid on
xlabel('x position (m)')
ylabel('y position (m)')
title('Trajectory of Electrons in Random Motion with Scaterring Probability')


figure(3)
plot(temp)
grid on
xlabel('Time')
ylabel('Temp (K)')
title('Semiconductor Temperature vs. Time')