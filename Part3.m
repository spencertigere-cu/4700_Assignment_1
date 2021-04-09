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
tmn = 2e-13;                        
freepath = thermal_vel*2e-13;

num_E = 50;
E_xpos = rand(1,num_E).*r_w;
E_ypos = rand(1,num_E).*r_l;

% I lost marks in my first submission because electrons were leaking into
% the box. This resolves the issue
outbox = 1;
for k = 1:num_E
        
         while (E_xpos(k) >= 0.8e-7 & E_xpos(k) <= 1.2e-7) & ((E_ypos(k) >= 0 & E_ypos(k) <= 0.4e-7)|(E_ypos(k) >= 0.6e-7 & E_ypos(k) <= 1e-7))
              E_xpos(k) = rand(1, 1).* r_w;
              E_ypos(k) = rand(1, 1).* r_l;
         end

end

angle = rand(1,num_E).*2*pi;

sig = sqrt(boltzmann_constant*T/effective_mass)/4;
MBdist = makedist('Normal',thermal_vel,sig);
E_vel = random(MBdist,1,num_E);
angle = rand(1,num_E).*2*pi;
E_xvel = E_vel.*cos(angle);
E_yvel = E_vel.*sin(angle);
d_t = 1e-9/thermal_vel;
p_scat = 1 - exp(-d_t/tmn);

for i=1:1000
      for j = 1:num_E 
          if p_scat > rand
                angleNew = rand(1).*2*pi;
                E_xvel(1,j) = random(MBdist,1).*cos(angleNew);
                E_yvel(1,j) = random(MBdist,1).*sin(angleNew);
          end
             if E_ypos(1,j) + E_yvel(1,j).*d_t >=1e-7|| E_ypos(1,j) + E_yvel(1,j).*d_t <= 0||((E_ypos(1,j) + E_yvel(1,j).*d_t >= 0.6e-7||E_ypos(1,j) + E_yvel(1,j).*d_t <= 0.4e-7) & (E_xpos(1,j) + E_xvel(1,j)*d_t>=0.8e-7 & E_xpos(1,j) + E_xvel(1,j)*d_t<=1.2e-7))
                E_yvel(1,j) = E_yvel(1,j)*-1;          
             end
             ypos(i,j) = E_ypos(1,j) + E_yvel(1,j).*d_t;
             if E_xpos(1,j) + E_xvel(1,j)*d_t >= r_w || E_xpos(1,j) + E_xvel(1,j)*d_t <= 0||((E_ypos(1,j) + E_yvel(1,j).*d_t >= 0.6e-7||E_ypos(1,j) + E_yvel(1,j).*d_t <= 0.4e-7) & (E_xpos(1,j) + E_xvel(1,j)*d_t>=0.8e-7 & E_xpos(1,j) + E_xvel(1,j)*d_t<=1.2e-7))
                if E_xpos(1,j) + E_xvel(1,j)*d_t >= r_w ||  E_xpos(1,j) + E_xvel(1,j)*d_t >= 0.8e-7
                    xpos(i,j) = 0;
                else
                    xpos(i,j) = r_w;
                end
             else
                xpos(i,j) = E_xpos(1,j) + E_xvel(1,j).*d_t;
             end
              t1 = sqrt(E_xvel(1,j).^2 + E_yvel(1,j).^2);
              temp(i,j) = ((mean(t1)^2)*effective_mass)/(2*boltzmann_constant);
              tempx(i,j) = E_xpos(1,j);
              tempy(i,j) = E_ypos(1,j);
        end
    
    E_xpos = xpos(i,:);
    E_ypos = ypos(i,:);
end

plot(xpos,ypos,'.','MarkerSize',12);
hold on
topboxy = [r_l 0.6e-7 0.6e-7 r_l];
boxx = [0.8e-7 0.8e-7 1.2e-7 1.2e-7];
bottomboxy = [0 0.4e-7 0.4e-7 0];
plot(boxx,topboxy,'LineWidth',4);
plot(boxx,bottomboxy,'LineWidth',4);
xlim([-0.1e-7 2.1e-7]);
ylim([-0.1e-7 1.1e-7]);
grid on
xlabel('x position (m)')
ylabel('y position (m)')
title('Trajectory of Electrons in Random Motion with Scaterring Probability')
hold off

figure(2)
E_den = [xpos(1000,:)',ypos(1000,:)'];

hist3(E_den,'CdataMode','auto')
xlabel('x position(m)')
ylabel('y position(m)')
title('Electron Density Map')
colorbar
view(2)

%I lost marks in my first submission because the histogram was missing.
%This has been fixed.
figure(3)
[X,Y] = meshgrid(0e-7:0.1e-7:2e-7,0e-7:0.1e-7:1e-7);
tempNew = griddata(tempx,tempy,temp,X,Y);
surf(X,Y,tempNew)
colorbar
xlabel('Width')
ylabel('Length')
title('Temperature Mapping')