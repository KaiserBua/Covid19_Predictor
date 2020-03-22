clear all
close all

data = readtable('JohnHopkins_21.03_Deutschland.csv');
neuInfizierte = data.difference;
neuInfizierte(neuInfizierte<0) = 0;
%neuInfizierte = neuInfizierte(34:end)

neuInfizierte = [2;3;5;27;13; 51;33;38;52;160;239;156;107;237;157;271;802;693;733;1043;1174;1144;1042;2801;2958;2705];
%plot(neuInfizierte)


for rueckSchau_Zeitraum = 1:20

%%anstiegFaelleInGerProTagHistory = [2;3;5;27;13; 51;33;38;52;160;239;156;107;237;157;271;802;693;733;1043;1174;1144;1042;2801;2958;2705];
cumulativeCases = zeros(length(neuInfizierte),1);
cumulativeCases(1) = neuInfizierte(1);
for i = 2:length(neuInfizierte)
    cumulativeCases(i) = cumulativeCases(i) + neuInfizierte(i);
end

indizes = length(neuInfizierte)-rueckSchau_Zeitraum + 1 : length(neuInfizierte); 

x= transpose(indizes);
y =log(cumulativeCases(indizes));
m(rueckSchau_Zeitraum) = x\y; 
verdopplungzeit = log(2)./m

end

figure(1) 
plot(verdopplungzeit) 

% figure(2)
% 
% subplot(1,2,1)
% hold on
% 
% plot(log(exp(indizes)))
% plot(log(cumulativeCases(indizes)));
% plot(m*x)
% axis equal
% 
% 
% subplot(1,2,2)
% hold on
% plot(cumulativeCases)
% x1 = length(neuInfizierte)-rueckSchau_Zeitraum + 1 : length(neuInfizierte); 
% plot(x1, cumulativeCases(12)*exp(m*(x1-12)));
% verdopplungzeit = log(2)/m