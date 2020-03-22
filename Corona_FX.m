close all
clear all

AnteilImmuner = 0;
Basisreproduktionszahl = 1.8; % Source: RKI gibt an, wie viele Menschen eine bereits erkrankte Person im Durchschnitt infiziert, falls die betroffene Bevölkerung weder geimpft noch anderweitig vor der Übertragung geschützt wird 
Nettoreproduktionszahl = Basisreproduktionszahl * (1-AnteilImmuner); % Source WIKI
Vorhersage = 150;

AnteilKrankerZuBehandelnder = 0.1;
AnteilBehandelterAufIntensivMitBeatmung = 0.1;
AnteilToterUnterSchwerenFaellen = 0.22;

Population = 80e6 * ones(1,Vorhersage);

Inkubationszeit = 5;
InfektioeseZeit = 4;
PneunmonieZeit = 3;
KrankenhausPhase = 14;

neuInfizierteAmTag = zeros(1,Vorhersage+8);
bisherInfizierteAmTag = zeros(1,Vorhersage);
geneseneAmTag = zeros(1,Vorhersage);
InkubierendeInfizierteAmTag = zeros(1, Vorhersage);
InfektioeseInfizierteAmTag = zeros(1,Vorhersage);
ZubehandelndeInfizierteAmTag = zeros(1,Vorhersage);
IntensivPatientenAmTag = zeros(1,Vorhersage);
ToteAmTag= zeros(1,Vorhersage);
bisherToteAmTag= zeros(1,Vorhersage);
neuGeneseneAmTag = zeros(1,Vorhersage);

neuInfizierteAmTag(1) = 10000;

for t = 1:Vorhersage
    
    % Gleichverteilte Infektion neuer Personen über den Infektionszeitraum
    for x = (t+Inkubationszeit):(t+Inkubationszeit+InfektioeseZeit-1)
       
        neuInfizierteAmTag(x) = neuInfizierteAmTag(x) + neuInfizierteAmTag(t)*Nettoreproduktionszahl / InfektioeseZeit;
       
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
                   
    
    
    AnteilImmuner = bisherInfizierteAmTag(t)/Population(t);
    
    Nettoreproduktionszahl = Basisreproduktionszahl * (1-AnteilImmuner);
    
    
end

figure(1)
hold on
plot(neuInfizierteAmTag,'r');
plot(bisherInfizierteAmTag,'b');
plot(geneseneAmTag, 'g');
plot(InkubierendeInfizierteAmTag,'y');
plot(Population-geneseneAmTag,'c');
legend('NeuInfizierteAmTag', 'bisherinfizierteAmTag','geneseneAmTag','InkubierendeInfizierteAmTag','Infizierbare Bevölkerung');
title('Infektionsgeschehen in der Bevölkerung');
figure(2)
hold on
plot(ZubehandelndeInfizierteAmTag,'g')
plot(IntensivPatientenAmTag,'y')
plot(ToteAmTag,'r')
plot(bisherToteAmTag,'bl')
legend('ZubehandelndeInfizierteAmTag','IntensivPatientenAmTag','ToteAmTag','bisherToteAmTag')
title('Auslastung Gesundheitssystem');