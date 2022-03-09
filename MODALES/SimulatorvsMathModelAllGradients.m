%HouseKeeping
clear;
clc;
close all;
%% Questions
% theta is in degrees
%%  Vehicle model variables
% Load Excel file data
load('thetadistance.mat');
Mv = 1505;%kg %Me
g = 9.81/1000; % km/s^2
Cr = 0.012;
Cd = 0.31;

% calculate average of every 20 values
grade = computeAverage(theta, 20); 
dd = computeAverage(distance, 20)./1000;  %ones(length(grade),1)*0.1;

rho = 1.2;
Af = 2/1000; % m^2
%Va = 18.75; % average speed in km/s (taken from Jianbing et al. paper)
dv = 0;
dt = 0;

%% Computing Feng and FC (using quadratic equation from Feng (power));
%adjmat = csvread('AdjacencyMatrix.csv');
adjmat = csvread('adjmatv2.csv');
temp = adjmat;
for i=1:length(grade)
    adjmat(:,:,i)=temp;
end
%dtt = load("averageTime.mat");
%dtt = dtt.dt;
%dtt(find(dtt==0))=Inf;  % make all zero distances equal to infinity
startSpeed= (0:10:130)'./3600; %km/s
endSpeed = (0:10:130)'./3600;%km/s
Feng = zeros(14,14, length(grade));
FC = zeros(14,14, length(grade));
NOX = zeros(14,14, length(grade));
power = zeros(14,14, length(grade));
dt = zeros(14,14, length(grade));
Fbrk = 0;
for index=1: length(grade)
    for i = 1:14
        for j = 1:14
            dv = (endSpeed(j,1)-startSpeed(i,1)); %km/s
            dt(i,j,index) =  2*(dd(index)./(startSpeed(i,1)+endSpeed(j,1)));  %s       %dtt(i,j);
            Va = (startSpeed(i,1)+endSpeed(j,1))/2; % km/s       %(dd(index)./dt); 
            Frol = (Mv*g*Cr*cosd(grade(index))); %kg km/s^2
            Faro = ((1/2)* rho* Af * Cd * Va^2)*1000;% 
            Fgrd = (Mv* g* sind(grade(index))); %kg km/s^2
            Ftotal(i,j,index) = Frol + Faro + Fgrd;

            Feng(i,j, index) = Mv*(dv./dt(i,j,index))+Ftotal(i,j,index); %kg km/s^2
            
            power(i,j, index) = Va.*Feng(i,j, index).*1000;% 1k Watt=kg m^2/s^3 https://en.wikipedia.org/wiki/Watt
            % NOx emmission costs
            %NOX(i,j,index) = ((0.0203*power(i,j, index)^2 + 0.2062*power(i,j, index) + 0.4204)); %milligrams/s
            %NOXX(i,j,index) = (NOX(i,j,index)./Va).*10^-3; %g/km
            %NOXX(find(NOXX<0))=0.4204;  % make all negative costs equal to 0.4204
            % NOx emmission costs
            if power(i,j,index)<-2 %kW
                NOX(i,j,index) = 0;
            else
                NOX(i,j,index) = ((0.0203*power(i,j, index)^2 + 0.2062*power(i,j, index) + 0.4204)); %milligrams/s
            end
            NOXX(i,j,index) = (NOX(i,j,index)./Va).*10^-3; %g/km
            % Fuel consumption costs
             if (power(i,j, index)>=0)
                %FC(i,j,index) =((2.595002839079077*10^-07)*Feng(i,j)^2) + ((5.400776287997336*10^-05)*Feng(i,j))+(1.322524181578653*10^-04);
                FC(i,j,index) = ((0.0009*power(i,j, index)^2)+(0.1944*power(i,j, index))+(0.4761));%kg/h
                FCC(i,j,index) = (FC(i,j,index)./3.6)./(Va);  %g/km
        
             else 
                 %FC(i,j,index)=(1.322524181578653*10^-04);
                 FC(i,j,index)= 0.4761;%kg/h
                 FCC(i,j,index) = (FC(i,j,index)./3.6)./(Va);%grams/meter
             end
        end
    end
end
   a =  498.4;  %(164.2, 832.7)
   b =  229.2;  %(-53.87, 512.2)
   c =  -0.5289;  %(-0.7702, -0.2876)
   power = power.*adjmat;
   Feng =  Feng.*adjmat;
   FC = FC.*adjmat;
   FCC = FCC.*adjmat;
   %FCCq = polyval(p1,FCC);
   FCCe = a.*exp(-(FCC./b).^c);
   FCC=FCCe;
   NOX = NOX.*adjmat;
   NOXX = NOXX.*adjmat;
   dt = dt.*adjmat;
%    [r, c] = size(FC);
    fuelcosts = csvread('fuelconsumption.csv');
    fccgrade0=load('Fuelcosts_grade0.mat');
    fccgrade0 = fccgrade0.FCC;
%    fccgrade3=load('FCCgrade3.mat');
%    fccgrade3 = fccgrade3.FCCgrade3;
   sz=120;
   figure(1);
   hold on;
   for i=1:14
       for j=1:14
           if(FCC(i,j)~=0)
            scatter(FCC(i,j,1),(i-1)*10,sz,'g+');%+ve gradient
            scatter(FCC(i,j,132),(i-1)*10,sz,'r*');%-ve gradient
            %scatter(fccgrade0(i,j),(i-1)*10,sz,'bo');
            scatter(fuelcosts(i,j),(i-1)*10,sz,'bo');

            %hold on;
           end
       end
   end
   legend('+ve gradient','-ve gradient','Grade 0');
   %legend('Simulator', 'Mathematical Model');
   ylabel('Start Speed (km/h)');
   xlabel('Costs (g/km)');
   set(gca,'XAxisLocation','top','YAxisLocation','left','ydir','reverse');
   %zlabel('Cost (g/km)');
   %ylim([0 1])
   hold off;
   title('Fuel Costs (g/km)');