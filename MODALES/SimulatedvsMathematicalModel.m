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
grade = 0; % degree %mean(theta(1:20,1));%computeAverage(theta, 20); 
dd =  0.1; % kilometers %mean(distance(1:20,1));%computeAverage(distance, 20); 

rho = 1.2;
Af = 2; % m^2
%Va = 18.75; % average speed in km/s (taken from Jianbing et al. paper)
dv = 0;
dt = 0;

%% Computing Feng and FC (using quadratic equation from Feng (power));
%adjmat = csvread('AdjacencyMatrix.csv');
adjmat = csvread('adjmatv2.csv');
adjmat(8,12)=1;
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
            %NOXX(find(NOXX<0))=0.4204;  % make threshold costs equal to 0.4204

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

%% Hierarchical Modeling %%
   %p1 = [-0.00243395450802555 1.44299084392367 0.306143141740943];
   p1 = [-0.0017 1.1931 -0.4611];
   p2 = [-0.1042 1.4861 0.0582];
   a =  498.4;  %(164.2, 832.7)
   b =  229.2;  %(-53.87, 512.2)
   c =  -0.5289;  %(-0.7702, -0.2876)
   an = 512.7;  %(-1.209e+04, 1.311e+04)
   bn = 5.9e+07;  %(-6.589e+09, 6.707e+09)
   cn = -0.1033;  %(-0.537, 0.3304)
   power = power.*adjmat;
   Feng =  Feng.*adjmat;
   FC = FC.*adjmat;
   FCC = FCC.*adjmat;
   FCCq = polyval(p1,FCC);
   FCCe = a.*exp(-(FCC./b).^c);
   FCC=FCCq;
   NOX = NOX.*adjmat;
   NOXX = NOXX.*adjmat;
   NOXXq = polyval(p2,NOXX);
   NOXXq = NOXXq.*adjmat;
   NOXXe = an.*exp(-(NOXX./bn).^cn);
   %NOXX = NOXXq;
   dt = dt.*adjmat;
   [r, c] = size(FC);
   fuelcosts = csvread('fuelconsumption.csv');
%    fccgrade0=load('Fuelcosts_grade0.mat');
%    fccgrade0 = fccgrade0.FCC;
%    fccgrade3=load('FCCgrade3.mat');
%    fccgrade3 = fccgrade3.FCCgrade3;
    FCC(isnan(FCC))=0;
   sz=120;
   error=0;
   figure(1);
   hold on;
   for i=1:14
       for j=1:14
           if(FCC(i,j)~=0 && fuelcosts(i,j)~=0)
            error=error+abs(FCC(i,j)-fuelcosts(i,j));
            scatter(FCC(i,j),(i-1)*10,sz,'rh');
            %hold on;
            scatter(fuelcosts(i,j),(i-1)*10,sz,'bo');
            
           end
       end
       %pause;
   end
   %legend('Simulator', 'Mathematical Model');
   ylabel('Start Speed (km/h)');
   xlabel('Fuel Costs (g/km)');
   set(gca,'XAxisLocation','top','YAxisLocation','left','ydir','reverse');
   legend('Mathematical Model', 'Simulator');
   %zlabel('Cost (g/km)');
   %ylim([0 1])
   %hold off;
   %title('Fuel Costs (g/km) Error:',error);
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

   NOXX(isnan(NOXX))=0;
   figure(2);   
   hold on;
   noxcosts = csvread('nox.csv');
   
   for i=1:14
       for j=1:14
           if(NOXX(i,j)~=0 && NOXX(i,j)~=inf)
            scatter(NOXX(i,j),(i-1)*10,sz,'rh')
            scatter(noxcosts(i,j),(i-1)*10,sz,'go')
            %scatter3((i-1)*10,(j-1)*10,NOXX(i,j),'b*');
            %scatter3((i-1)*10,(j-1)*10,noxcosts(i,j),'r*');
            %noxdifference(i,j)=(NOXX(i,j)-noxcosts(i,j));
            %rnoxdifference(i,j)=(NOXX(i,j)-noxcosts(i,j))./NOXX(i,j);
            
           end
       end
   end
   ylabel('Start Speed (km/h)');
   xlabel('NOx Cost (g/km)');
   set(gca,'XAxisLocation','top','YAxisLocation','left','ydir','reverse');
   legend('Mathematical Model','Simulator');
   %xlim([0 6])
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