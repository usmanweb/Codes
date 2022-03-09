%housekeeping
clc;
clear;
close all;
a=0:0.001:200;
b=0:0.001:200;
sz=120;
noxcosts = csvread('nox.csv');%Simulator
noxgrade0=load('Noxx.csv');%MathModel
figure(1);
hold on;
noxgrade0u=triu(noxgrade0);
noxcostsu=triu(noxcosts);
for i=1:length(noxgrade0u)
    for j=1:length(noxgrade0u)
        
        if noxgrade0u(i, j) ~= 0 && noxcostsu(i, j) ~= 0
            plot(noxgrade0u(i,j),noxcostsu(i,j),'r*','LineWidth',2);
        else
            %plot(noxgrade0u(i,j),noxcostsu(i,j),'b*','LineWidth',2);
        end
    end
end
ylabel('Simulator nox costs');
xlabel('Mathematical Model nox costs');
title('Upper diagonal matrix (Speeding up)');

figure(2);
hold on;
noxgrade0l=tril(noxgrade0,-1);
noxcostsl=tril(noxcosts,-1);
for i=1:length(noxgrade0l)
    for j=1:length(noxgrade0l)
        
        if noxgrade0l(i, j) ~= 0 && noxcostsl(i, j) ~= 0
            plot(noxgrade0l(i,j),noxcostsl(i,j),'b*','LineWidth',2);
        else
            %plot(noxgrade0u(i,j),noxcostsu(i,j),'b*','LineWidth',2);
        end
    end
end
ylabel('Simulator nox costs');
xlabel('Mathematical Model nox costs');
title('Lower diagonal matrix (Slowing down)')