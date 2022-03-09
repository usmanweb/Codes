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
grade = -3; % degree %mean(theta(1:20,1));%computeAverage(theta, 20); 
dd =  0.1; % kilometers %mean(distance(1:20,1));%computeAverage(distance, 20); 

rho = 1.2;
Af = 2/1000; % m^2
%Va = 18.75; % average speed in km/s (taken from Jianbing et al. paper)
dv = 0;
dt = 0;

%% Computing Feng and FC (using quadratic equation from Feng (power));
%adjmat = csvread('AdjacencyMatrix.csv');
adjmat = csvread('adjmatv2.csv');
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
            Ftotal = Frol + Faro + Fgrd;

            Feng(i,j, index) = Mv*(dv./dt(i,j,index))+Ftotal; %kg km/s^2
            
            power(i,j, index) = Va.*Feng(i,j, index).*1000;% 1k Watt=kg m^2/s^3 https://en.wikipedia.org/wiki/Watt
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
   [r, c] = size(FC);
   fuelcosts = csvread('fuelconsumption.csv');
   fccgrade0=load('Fuelcosts_grade0.mat');
   fccgrade0 = fccgrade0.FCC;
   fccgrade3=load('FCCgrade3.mat');
   fccgrade3 = fccgrade3.FCCgrade3;
   sz=120;
   figure(1);
   for i=1:14
       for j=1:14
           if(FCC(i,j)~=0)
            scatter(fccgrade3(i,j),(i-1)*10,sz,'g+');
            scatter(FCC(i,j),(i-1)*10,sz,'r*');
            scatter(fccgrade0(i,j),(i-1)*10,sz,'bo');

%             scatter3((i-1)*10,(j-1)*10,fccgrade3(i,j),sz,'g+');
%             scatter3((i-1)*10,(j-1)*10,FCC(i,j),sz,'r*');
%             scatter3((i-1)*10,(j-1)*10,fccgrade0(i,j),sz,'bo');
            %scatter(fuelcosts(i,j),(i-1)*10,sz,'bo');
            %scatter3((i-1)*10,(j-1)*10,FCC(i,j),'b*');
            %scatter3((i-1)*10,(j-1)*10,fuelcosts(i,j),'r*');
            %fueldifference(i,j)=(FCC(i,j)-fuelcosts(i,j));%./FCC(i,j);
            %rfueldifference(i,j)=(FCC(i,j)-fuelcosts(i,j))./FCC(i,j);
            %scatter((j-1)*10,fueldifference(i,j),'r*');
            hold on;
           end
       end
       %pause;
   end
   legend('Grade 0', 'Grade +3','Grade -3');
   %legend('Simulator', 'Mathematical Model');
   ylabel('Start Speed (km/h)');
   xlabel('Costs (g/km)');
   set(gca,'XAxisLocation','top','YAxisLocation','left','ydir','reverse');
   %zlabel('Cost (g/km)');
   %ylim([0 1])
   hold off;
   title('Fuel Costs (g/km)');
%    figure(2);
%    Z = fueldifference;
%    %surf(Z)
%    %hold on
%    imAlpha=ones(size(Z));
%    imAlpha(isnan(Z))=0;
%    imAlpha(find(Z==0))=0;
%    imAlpha(find(Z==1))=0;
%    imagesc(Z,'AlphaData',imAlpha);
%    colorbar;
%    title('Absolute Error Difference for Fuel values');
%    figure(3);
%    imAlpha=ones(size(Z));
%    imAlpha(isnan(Z))=0;
%    imAlpha(find(Z==0))=0;
%    imAlpha(find(Z==1))=0;
%    imagesc(rfueldifference,'AlphaData',imAlpha);
%    colorbar;
%    title('Relative Error Difference for Fuel values');

   
%    figure(4);   
%     
%    noxcosts = csvread('nox.csv');
%    
%    for i=1:14
%        for j=1:14
%            if(NOXX(i,j)~=0 && NOXX(i,j)~=inf)
%             scatter((i-1)*10,NOXX(i,j),'b*')
%             scatter((i-1)*10,noxcosts(i,j),'d')
%             %scatter3((i-1)*10,(j-1)*10,NOXX(i,j),'b*');
%             %scatter3((i-1)*10,(j-1)*10,noxcosts(i,j),'r*');
%             %noxdifference(i,j)=(NOXX(i,j)-noxcosts(i,j));
%             %rnoxdifference(i,j)=(NOXX(i,j)-noxcosts(i,j))./NOXX(i,j);
%             hold on;
%            end
%        end
%    end
%    legend('Simulator', 'Mathematical Model');
%    xlabel('Start Speed (km/h)');
%    ylabel('End Speed(km/h)');
%    zlabel('Cost (g/km)');
%    %ylim([0 1])
%    hold off;
%    title('NOx Costs (g/km)');
%    figure(5);
%    Znox = noxdifference;
%    %surf(Z)
%    %hold on
%    imAlpha=ones(size(Znox));
%    imAlpha(isnan(Znox))=0;
%    imAlpha(find(Znox==0))=0;
%    imAlpha(find(Znox==1))=0;
%    imagesc(Znox,'AlphaData',imAlpha);
%    colorbar;
%    figure(6);
%    Znoxx = rnoxdifference;
%    %surf(Z)
%    %hold on
%    imAlpha=ones(size(Znoxx));
%    imAlpha(isnan(Znoxx))=0;
%    imAlpha(find(Znoxx==0))=0;
%    imAlpha(find(Znoxx==1))=0;
%    imagesc(Znoxx,'AlphaData',imAlpha);
%    colorbar;