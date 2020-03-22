close all
clear all


%%Inspired by Christopher Roemmelmayer
%%Implemented by Felix Kaiser

AnteilImmuner = 0;


Basisreproduktionszahl = 1.8; % Source: RKI gibt an, wie viele Menschen eine bereits erkrankte Person im Durchschnitt infiziert, falls die betroffene Bevölkerung weder geimpft noch anderweitig vor der Übertragung geschützt wird 
Nettoreproduktionszahl = Basisreproduktionszahl * (1-AnteilImmuner); % Source WIKI
Vorhersagezeitraum = 300;

AnteilKrankerZuBehandelnder = 0.1;% Konservative Wahrscheinlichkeiten, behandelt werden zu müssen
AnteilBehandelterAufIntensivMitBeatmung = 0.1;%Intensiv behandelt werden zu müssen
AnteilToterUnterSchwerenFaellen = 0.22; % zu sterben



Inkubationszeit = 5; %Source RKI Steckbrief
InfektioeseZeit = 4;
PneunmonieZeit = 3;
KrankenhausPhase = 14;


Population = 80e6 * ones(1,Vorhersagezeitraum);
neuInfizierteAmTag = zeros(1,Vorhersagezeitraum+8);
bisherInfizierteAmTag = zeros(1,Vorhersagezeitraum);
geneseneAmTag = zeros(1,Vorhersagezeitraum);
InkubierendeInfizierteAmTag = zeros(1, Vorhersagezeitraum);
InfektioeseInfizierteAmTag = zeros(1,Vorhersagezeitraum);
ZubehandelndeInfizierteAmTag = zeros(1,Vorhersagezeitraum);
IntensivPatientenAmTag = zeros(1,Vorhersagezeitraum);
ToteAmTag= zeros(1,Vorhersagezeitraum);
bisherToteAmTag= zeros(1,Vorhersagezeitraum);
neuGeneseneAmTag = zeros(1,Vorhersagezeitraum);
AnteilImmuner = zeros(1,Vorhersagezeitraum);



anstiegFaelleInGerProTagHistory = [2;3;5;27;13; 51;33;38;52;160;239;156;107;237;157;271;802;693;733;1043;1174;1144;1042;2801;2958];

simStart = length(anstiegFaelleInGerProTagHistory);

neuInfizierteAmTag(1:simStart) = anstiegFaelleInGerProTagHistory;

est_Nettoreproduktionszahl = anstiegFaelleInGerProTagHistory(end)*InfektioeseZeit / sum(anstiegFaelleInGerProTagHistory(end-(Inkubationszeit+InfektioeseZeit):end-Inkubationszeit))

BasisreproduktionsrateArray = ones(1,Vorhersagezeitraum);
NettoreproduktionszahlArray = ones(1,Vorhersagezeitraum);

%%%%CHANGE MODE AND PARAMETERS HERE
MODE = 3; 

if (MODE == 1)
    disp('Manuelle Basisreproduktionszahl');
    BasisreproduktionsrateArray = ones(1,Vorhersagezeitraum) * Basisreproduktionszahl;
    NettoreproduktionszahlArray = ones(1,Vorhersagezeitraum) * Basisreproduktionszahl;
elseif (MODE == 2)
    
    disp('Berechnete Basisreproduktionszahl');
    est_Nettoreproduktionszahl = anstiegFaelleInGerProTagHistory(end)*InfektioeseZeit / sum(anstiegFaelleInGerProTagHistory(end-(Inkubationszeit+InfektioeseZeit):end-Inkubationszeit))
    disp(est_Nettoreproduktionszahl);
    BasisreproduktionsrateArray = ones(1,Vorhersagezeitraum) * est_Nettoreproduktionszahl;
    NettoreproduktionszahlArray = ones(1,Vorhersagezeitraum) * est_Nettoreproduktionszahl;
        
elseif (MODE == 3)
        
    disp('Manuelle Basisreproduktionszahl + Ausgangssperre mit einstellbarer Ausgangssperre');
    
    Tag_der_Ausgangssperre = 30;
    Tag_des_Ende_der_Ausgangssperre = 51;
    Basisreproduktionsrate_Nach_Ausgangssperre = 0.8;
        
    BasisreproduktionsrateArray = ones(1,Vorhersagezeitraum);
    BasisreproduktionsrateArray(1:end) = Basisreproduktionszahl;
    BasisreproduktionsrateArray(Tag_der_Ausgangssperre + 1:Tag_des_Ende_der_Ausgangssperre) = Basisreproduktionsrate_Nach_Ausgangssperre;
    
    NettoreproduktionszahlArray =  ones(1,Vorhersagezeitraum);
end


IMPFUNG = 0;

if (IMPFUNG)
    
    disp('Impfung aktiviert');
    Tag_der_Markteinfuehrung = 50;
    taegliche_Impfungen = 5e5;
end


GesamtIntensivBetten = 28000; % Gesamt Intensivbetten
IntensivKapazitaet = zeros(1,Vorhersagezeitraum) + GesamtIntensivBetten*0.21;

INTENSIV_AUSBAU = 0;

if (INTENSIV_AUSBAU)
    Tag_des_Beginns = simStart;
    taeglich_neue_Betten = 1000;
    
    for ausbauTag = (Tag_des_Beginns:Vorhersagezeitraum)
        IntensivKapazitaet(ausbauTag) = IntensivKapazitaet(ausbauTag-1) + taeglich_neue_Betten;

    end
    
end











for t = 1:Vorhersagezeitraum
    
    if t>=simStart-Inkubationszeit
        % Gleichverteilte Infektion neuer Personen über den Infektionszeitraum
        for x = (t+Inkubationszeit):(min(Vorhersagezeitraum,t+Inkubationszeit+InfektioeseZeit-1))

            neuInfizierteAmTag(x) = neuInfizierteAmTag(x) + neuInfizierteAmTag(t)*NettoreproduktionszahlArray(t) / InfektioeseZeit;

        end
    end
    
    %Infizierte in Inkubationszeit
    InkubierendeInfizierteAmTag(t) = sum(neuInfizierteAmTag(max(1,t-(Inkubationszeit)):t));
                         
    %Infektiöse und symptomezeigende Infizierte
    InfektioeseInfizierteAmTag(t) = sum(neuInfizierteAmTag(max(1,t-(InfektioeseZeit+Inkubationszeit)):t-Inkubationszeit));
    
    %ZubehandelndeInfizierte
    ZubehandelndeInfizierteAmTag(t) = AnteilKrankerZuBehandelnder * sum(neuInfizierteAmTag(max(1,t-(Inkubationszeit+InfektioeseZeit+PneunmonieZeit)):max(1,t-(Inkubationszeit+InfektioeseZeit))));
  
    %IntensivPatienten
    IntensivPatientenAmTag(t) = AnteilKrankerZuBehandelnder * AnteilBehandelterAufIntensivMitBeatmung*  sum(neuInfizierteAmTag(max(1,t-(Inkubationszeit+InfektioeseZeit+PneunmonieZeit+KrankenhausPhase)):max(1,t-(Inkubationszeit+InfektioeseZeit+PneunmonieZeit))));
    
    %Tote
    if t>Inkubationszeit+InfektioeseZeit+PneunmonieZeit+KrankenhausPhase
        ToteAmTag(t) = AnteilKrankerZuBehandelnder * AnteilBehandelterAufIntensivMitBeatmung * AnteilToterUnterSchwerenFaellen * neuInfizierteAmTag(t-(Inkubationszeit+InfektioeseZeit+PneunmonieZeit+KrankenhausPhase));
    end
    bisherToteAmTag(t) = bisherToteAmTag(max(1,t-1)) + ToteAmTag(t);
    Population(t) = Population(max(1,t-1))-ToteAmTag(t);
    
    %bisherInfizierte kumulativ
    bisherInfizierteAmTag(t) = bisherInfizierteAmTag(max(1,t-1)) + neuInfizierteAmTag(t);
    
    
    %Begrenzung: Es kann nur die ganze Population angesteckt werden
    if (bisherInfizierteAmTag(t) + neuInfizierteAmTag(t+1) > Population(t))
        neuInfizierteAmTag(t+1) = Population(t) - bisherInfizierteAmTag(t);
    end
    
    
    
    if(t>Inkubationszeit+InfektioeseZeit)
       % geneseneAmTag(t) = sum(neuInfizierteAmTag(1:t-(Inkubationszeit+InfektioeseZeit)));
    end
    
    if(t>Inkubationszeit+InfektioeseZeit)
        
        neuGeneseneAmTag(t) = (1-AnteilKrankerZuBehandelnder) * neuInfizierteAmTag(t-(Inkubationszeit+InfektioeseZeit));
        
        if(t>Inkubationszeit+InfektioeseZeit+ PneunmonieZeit)
            
            neuGeneseneAmTag(t) = neuGeneseneAmTag(t)+ AnteilKrankerZuBehandelnder *  (1- AnteilBehandelterAufIntensivMitBeatmung) * neuInfizierteAmTag(t-(Inkubationszeit+InfektioeseZeit + PneunmonieZeit));
    
    
            if (t>Inkubationszeit+InfektioeseZeit+ PneunmonieZeit+KrankenhausPhase)
                    
                 neuGeneseneAmTag(t) = neuGeneseneAmTag(t)+ AnteilKrankerZuBehandelnder *  AnteilBehandelterAufIntensivMitBeatmung * neuInfizierteAmTag(t-(Inkubationszeit+InfektioeseZeit + PneunmonieZeit + KrankenhausPhase)) - ToteAmTag(t);
    
            end    
        end
    end
    
    geneseneAmTag(t) = geneseneAmTag(max(1,t-1))+ neuGeneseneAmTag(t);
                   
    
    
    AnteilImmuner(t+1) = AnteilImmuner(t)+neuInfizierteAmTag(t)/Population(t);
    
    if (IMPFUNG == 1)
        if (t>Tag_der_Markteinfuehrung)
            AnteilImmuner(t+1) = AnteilImmuner(t+1) + taegliche_Impfungen/Population(t);
        end
    end
    
    if(AnteilImmuner(t+1) >1) 
        AnteilImmuner(t+1)=1;
    end
    
    NettoreproduktionszahlArray(t+1) = BasisreproduktionsrateArray(t) * (1-AnteilImmuner(t+1));
    
    
end

figure(1)
hold on
plot(neuInfizierteAmTag,'r');
plot(bisherInfizierteAmTag,'b');
plot(geneseneAmTag, 'g');
plot(InkubierendeInfizierteAmTag,'y');
plot(1:Vorhersagezeitraum,Population(1:Vorhersagezeitraum).*(1-AnteilImmuner(1:Vorhersagezeitraum)),'c');
title('Infektionsgeschehen in der Bevölkerung');
y = ylim; % current y-axis limits
plot([simStart simStart],[y(1) y(2)])
if(MODE == 3)
    plot([Tag_der_Ausgangssperre Tag_der_Ausgangssperre],[y(1) y(2)])
    plot([Tag_des_Ende_der_Ausgangssperre Tag_des_Ende_der_Ausgangssperre],[y(1) y(2)])
    legend('NeuInfizierteAmTag', 'bisherinfizierteAmTag','geneseneAmTag','InkubierendeInfizierteAmTag','Infizierbare Bevölkerung', 'Heute', 'Tag der Ausgangssperre', 'Ende der Ausgangssperre');
else
    legend('NeuInfizierteAmTag', 'bisherinfizierteAmTag','geneseneAmTag','InkubierendeInfizierteAmTag','Infizierbare Bevölkerung', 'Heute');
end
legend('Location','NW');
set(gcf,'Position',[10 10 2000 1000])
figure(2)
hold on
plot(ZubehandelndeInfizierteAmTag,'g')
plot(IntensivPatientenAmTag,'y')
plot(ToteAmTag,'r')
plot(bisherToteAmTag,'bl')

title('Auslastung Gesundheitssystem');
x = xlim; % current y-axis limits
plot(IntensivKapazitaet)


y = ylim; % current y-axis limits
plot([simStart simStart],[y(1) y(2)])

if(MODE == 3)
    plot([Tag_der_Ausgangssperre Tag_der_Ausgangssperre],[y(1) y(2)])
    plot([Tag_des_Ende_der_Ausgangssperre Tag_des_Ende_der_Ausgangssperre],[y(1) y(2)])
   legend('ZubehandelndeInfizierteAmTag','IntensivPatientenAmTag','ToteAmTag','bisherToteAmTag','Intensivbetten','Heute', 'Tag der Ausgangssperre', 'Ende der Ausgangssperre');
else
    legend('ZubehandelndeInfizierteAmTag','IntensivPatientenAmTag','ToteAmTag','bisherToteAmTag','Intensivbetten', 'Heute')
end
legend('Location','NW');
set(gcf,'Position',[10 10 2000 1000])
figure(3)
hold on
plot(neuInfizierteAmTag,'r');
plot(bisherInfizierteAmTag,'b');
plot(geneseneAmTag, 'g');
plot(InkubierendeInfizierteAmTag,'y');
plot(1:Vorhersagezeitraum,Population(1:Vorhersagezeitraum).*(1-AnteilImmuner(1:Vorhersagezeitraum)),'c');
title('Infektionsgeschehen in der Bevölkerung ZOOMED');
ylim([0 bisherInfizierteAmTag(end)])
y = ylim; % current y-axis limits
plot([simStart simStart],[y(1) y(2)])
if(MODE == 3)
    plot([Tag_der_Ausgangssperre Tag_der_Ausgangssperre],[y(1) y(2)])
    
    plot([Tag_des_Ende_der_Ausgangssperre Tag_des_Ende_der_Ausgangssperre],[y(1) y(2)])
    legend('NeuInfizierteAmTag', 'bisherinfizierteAmTag','geneseneAmTag','InkubierendeInfizierteAmTag','Infizierbare Bevölkerung', 'Heute', 'Tag der Ausgangssperre', 'Ende der Ausgangssperre');
else
    legend('NeuInfizierteAmTag', 'bisherinfizierteAmTag','geneseneAmTag','InkubierendeInfizierteAmTag','Infizierbare Bevölkerung', 'Heute');
end
legend('Location','NW');
set(gcf,'Position',[10 10 2000 1000])

figure(4)
hold on
plot(ZubehandelndeInfizierteAmTag,'g')
plot(IntensivPatientenAmTag,'y')
plot(ToteAmTag,'r')
plot(bisherToteAmTag,'bl')

title('Auslastung Gesundheitssystem ZOOMED');
xlim([0 Vorhersagezeitraum]);
ylim([0 1.1*IntensivKapazitaet(end)]);
plot(IntensivKapazitaet)


y = ylim; % current y-axis limits
plot([simStart simStart],[y(1) y(2)])

if(MODE == 3)
    plot([Tag_der_Ausgangssperre Tag_der_Ausgangssperre],[y(1) y(2)])
    
    plot([Tag_des_Ende_der_Ausgangssperre Tag_des_Ende_der_Ausgangssperre],[y(1) y(2)])
   legend('ZubehandelndeInfizierteAmTag','IntensivPatientenAmTag','ToteAmTag','bisherToteAmTag','Intensivbetten','Heute', 'Tag der Ausgangssperre', 'Ende der Ausgangssperre');
else
    legend('ZubehandelndeInfizierteAmTag','IntensivPatientenAmTag','ToteAmTag','bisherToteAmTag','Intensivbetten', 'Heute')
end
legend('Location','NW');
set(gcf,'Position',[10 10 2000 1000])
figure(5)
hold on
plot(neuInfizierteAmTag,'r');
plot(bisherInfizierteAmTag,'b');
plot(geneseneAmTag, 'g');
plot(InkubierendeInfizierteAmTag,'y');
plot(1:Vorhersagezeitraum,Population(1:Vorhersagezeitraum).*(1-AnteilImmuner(1:Vorhersagezeitraum)),'c');
title('Infektionsgeschehen in der Bevölkerung ZOOMED');
ylim([0 bisherInfizierteAmTag(simStart+10)])
y = ylim; % current y-axis limits
xlim([0 simStart+10]);
plot([simStart simStart],[y(1) y(2)])
if(MODE == 3)
    plot([Tag_der_Ausgangssperre Tag_der_Ausgangssperre],[y(1) y(2)])
    
    plot([Tag_des_Ende_der_Ausgangssperre Tag_des_Ende_der_Ausgangssperre],[y(1) y(2)])
    legend('NeuInfizierteAmTag', 'bisherinfizierteAmTag','geneseneAmTag','InkubierendeInfizierteAmTag','Infizierbare Bevölkerung', 'Heute', 'Tag der Ausgangssperre', 'Ende der Ausgangssperre');
else
    legend('NeuInfizierteAmTag', 'bisherinfizierteAmTag','geneseneAmTag','InkubierendeInfizierteAmTag','Infizierbare Bevölkerung', 'Heute');
end

legend('Location','NW');
set(gcf,'Position',[10 10 2000 1000])
