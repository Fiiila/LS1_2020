%% setup sequence
close all
clear all
clc

%% zadane parametry
Q1=1.5e-4;   %pritok [m^3/s]
H1=0.6;      %zadana pozadovana vyska hladiny v prvni nadrzi
H2=0.4;      %zadana pozadovana vyska hladiny v druhe nadrzi
S=25e-4;     %plocha dna jedne nadrze (plati pro obe)
cp=0.6;      %zadane konstanty
c2=cp;       %zadane konstanty
g=10;        %zaokrouhlene gravitacni zrychleni

%% a) zaroven vypocet parametru pro bod b)
%% vypocet nastaveni ventilu Sp a S2 pro zadany pritok

vp=sqrt(2*g*(H1-H2));   %rychlost proudeni kapaliny prepustecim ventilem
v2=sqrt(2*g*H2);        %rychlost proudeni kapaliny vytokem

Qp_1=Q1;      %podminka, pro ustalene hladiny je, ze prutok celym systemem se musi rovnat
Q2_1=Qp_1;

Sp_1=Qp_1/(cp*vp);  %obsah prepousteciho ventilu
S2_1=Q2_1/(c2*v2);  %obsah vyposteciho ventilu

%% vypocet nastaveni ventilu Sp a S2 pro zvetseny pritok

Qp_2=1.2*Q1;      %podminka, pro ustalene hladiny je, ze prutok celym systemem se musi rovnat
Q2_2=Qp_2;

Sp_2=Qp_2/(cp*vp);  %obsah prepousteciho ventilu
S2_2=Q2_2/(c2*v2);  %obsah vyposteciho ventilu

%% simulace nelinearniho modelu (prechodova charakteristika pro U=Q1)
%nastaveni podminek pro simulaci
t=250;          %delka (cas) simulace
dQ=0*Q1;        %velikost skoku na vstupu
stepTime=0;     %cas, ve kterem se udela skok o dQ
h1condition=0;  %pocatecni podminka pro H1
h2condition=0;  %pocatecni podminka pro H2
U=Q1;           %vstupni signal systemu
Sp=Sp_1;
S2=S2_1;
%pro pritok 1.2*Q1
%Sp=Sp_1;
%S2=S2_1;

%simulace
simoutNelinear=sim('dve_nadrze_nelinear');

%% Vykresleni simulace nelinearniho modelu (prechodova charakteristika pro U=Q1)
figure
plot(simoutNelinear.simout.Time,simoutNelinear.simout.Data(:,1))
title('Prechodova charakteristika nel. systemu pro pritok Q1')
xlabel('cas [s]')
ylabel('vyska hladiny [m]')
hold on
plot(simoutNelinear.simout.Time,simoutNelinear.simout.Data(:,2))
legend('prubeh H1','prubeh H2')

%% b)

%% Linearizace nelinearniho modelu pro zadany pritok Q1

Q_1=Q1;    %pritok pro reseni linearizace

%pro zadany pritok musim dopocist ustalene hodnoty H1r a H2r vyjadrenim z
%diff rovnic, kde je dH1 a dH2 =0
H2r_1=(Q_1/(c2*S2_1*sqrt(2*g)))^2; %ustaleny stav H2 pri pritoku Q
H1r_1=(Q_1/(cp*Sp_1*sqrt(2*g)))^2+(Q_1/(c2*S2_1*sqrt(2*g)))^2; %ustaleny stav H1 pri pritoku Q

%rovnice pro kontrolu...po dosazeni musi byt rovne nule (nebo co nejmensi)
diffH1_1=-(cp*Sp_1*sqrt(2*g*(H1r_1-H2r_1)))/S+Q_1/S;
diffH2_1=(cp*Sp_1*sqrt(2*g*(H1r_1-H2r_1)))/S-(c2*S2_1*sqrt(2*g*H2r_1))/S;

%maticeA
a11_1=-(cp*Sp_1*sqrt(2*g))/(2*S*sqrt(H1r_1-H2r_1));
a12_1=(cp*Sp_1*sqrt(2*g))/(2*S*sqrt(H1r_1-H2r_1));
a21_1=(cp*Sp_1*sqrt(2*g))/(2*S*sqrt(H1r_1-H2r_1));
a22_1=-(cp*Sp_1*sqrt(2*g))/(2*S*sqrt(H1r_1-H2r_1))-(c2*S2_1*sqrt(2*g))/(2*S*sqrt(H2r_1));
A_1=[a11_1,a12_1;a21_1,a22_1];

%matice B
B_1=[1/S;0];

%matice C
C_1=[1,0;0,1];

%matice D
D_1=[0;0];

%% Linearizace nelinearniho modelu pro ZVETSENY pritok 1.2*Q1

Q_2=1.2*Q1;    %pritok pro reseni linearizace

%pro zadany pritok musim dopocist ustalene hodnoty H1r a H2r vyjadrenim z
%diff rovnic, kde je dH1 a dH2 =0
H2r_2=(Q_2/(c2*S2_2*sqrt(2*g)))^2; %ustaleny stav H2 pri pritoku Q
H1r_2=(Q_2/(cp*Sp_2*sqrt(2*g)))^2+(Q_2/(c2*S2_2*sqrt(2*g)))^2; %ustaleny stav H1 pri pritoku Q

%rovnice pro kontrolu...po dosazeni musi byt rovne nule (nebo co nejmensi)
diffH1_2=-(cp*Sp_2*sqrt(2*g*(H1r_2-H2r_2)))/S+Q_2/S;
diffH2_2=(cp*Sp_2*sqrt(2*g*(H1r_2-H2r_2)))/S-(c2*S2_2*sqrt(2*g*H2r_2))/S;

%maticeA
a11_2=-(cp*Sp_2*sqrt(2*g))/(2*S*sqrt(H1r_2-H2r_2));
a12_2=(cp*Sp_2*sqrt(2*g))/(2*S*sqrt(H1r_2-H2r_2));
a21_2=(cp*Sp_2*sqrt(2*g))/(2*S*sqrt(H1r_2-H2r_2));
a22_2=-(cp*Sp_2*sqrt(2*g))/(2*S*sqrt(H1r_2-H2r_2))-(c2*S2_2*sqrt(2*g))/(2*S*sqrt(H2r_2));
A_2=[a11_2,a12_2;a21_2,a22_2];

%matice B
B_2=[1/S;0];

%matice C
C_2=[1,0;0,1];

%matice D
D_2=[0;0];

%% c

%% simulace - porovnani nelinearniho a linearizovanych modelu pro Q1 a 1.2*Q1

%nastaveni podminek pro simulaci nelinearni model
U=Q1;               %vstupni signal systemu
t=500;              %delka (cas) simulace
dQ=0.05*U;          %velikost skoku na vstupu
stepTime=t/4;       %cas, ve kterem se udela skok o dQ
h1condition=H1;     %pocatecni podminka pro H1
h2condition=H2;     %pocatecni podminka pro H2
H1r=h1condition;    %hladina ve ktere se system nachazi
H2r=h2condition;    %hladina ve ktere se system nachazi
%simulace nelinearniho modelu
nelinearSim=sim('dve_nadrze_nelinear');
%nastaveni parametru systemu 
A=A_1;
B=B_1;
C=C_1;
D=D_1;
U=0;
linearized_1Sim=sim('dve_nadrze_linearized');
A=A_2;
B=B_2;
C=C_2;
D=D_2;
linearized_2Sim=sim('dve_nadrze_linearized');

figure
subplot(2,1,1)
plot(nelinearSim.simout.Time,nelinearSim.simout.Data(:,1))
hold on
plot(linearized_1Sim.simout.Time,linearized_1Sim.simout.Data(:,1))
hold on
plot(linearized_2Sim.simout.Time,linearized_2Sim.simout.Data(:,1))
title('Odezva na skok dQ=0.05*Q1 v case t/4')
xlabel('cas [s]')
ylabel('H1 [m]')
grid on
legend('nelinearni','linearizovany pro Q1','linearizovany pro 1.2*Q1','Location','southeast')

subplot(2,1,2)
plot(nelinearSim.simout.Time,nelinearSim.simout.Data(:,2))
hold on
plot(linearized_1Sim.simout.Time,linearized_1Sim.simout.Data(:,2))
hold on
plot(linearized_2Sim.simout.Time,linearized_2Sim.simout.Data(:,2))
title('Odezva na skok dQ=0.05*Q1 v case t/4')
xlabel('cas [s]')
ylabel('H2 [m]')
grid on
legend('nelinearni','linearizovany pro Q1','linearizovany pro 1.2*Q1','Location','southeast')

%% staticka charakteristika
length=100;%pocet vzorku..kolikrat ma probehnout jednotlive simulace
Uhodnoty=linspace(0,2*Q1,length);
%inicializace pole pro vysledky simulaci
statNelinearH1=linspace(0,0,length);
statLin_1H1=linspace(0,0,length);
statLin_2H1=linspace(0,0,length);
statNelinearH2=linspace(0,0,length);
statLin_1H2=linspace(0,0,length);
statLin_2H2=linspace(0,0,length);
pocet=1;
%cyklus vypocitavajici potrebna data pro statickou charaktteristiku
for i = Uhodnoty
    %nastaveni podminek pro simulaci nelinearni model
    U=i;               %vstupni signal systemu
    t=2000;              %delka (cas) simulace
    dQ=0*U;          %velikost skoku na vstupu
    stepTime=t/4;       %cas, ve kterem se udela skok o dQ
    h1condition=0;     %pocatecni podminka pro H1
    h2condition=0;     %pocatecni podminka pro H2
    H1r=h1condition;    %hladina ve ktere se system nachazi
    H2r=h2condition;    %hladina ve ktere se system nachazi
    %simulace nelinearniho modelu
    statNelinear=sim('dve_nadrze_nelinear');
    %nastaveni parametru systemu 
    A=A_1;
    B=B_1;
    C=C_1;
    D=D_1;
    statLin_1=sim('dve_nadrze_linearized');
    U=1.2*i;%zvyseni pritoku pro model linearizovany pro zvyseny pritok
    A=A_2;
    B=B_2;
    C=C_2;
    D=D_2;
    statLin_2=sim('dve_nadrze_linearized');
    %zapiseme posledni (ustalene) hodnoty do pripravenych poli (predpokladame, ze delka simulace je dostatecna, aby se vystupy ustalily)
    statNelinearH1(pocet)=statNelinear.simout.Data(end,1);
    statLin_1H1(pocet)=statLin_1.simout.Data(end,1);
    statLin_2H1(pocet)=statLin_2.simout.Data(end,1);
    statNelinearH2(pocet)=statNelinear.simout.Data(end,2);
    statLin_1H2(pocet)=statLin_1.simout.Data(end,2);
    statLin_2H2(pocet)=statLin_2.simout.Data(end,2);
    pocet=pocet+1;
end
%% Vykresleni staticke charakteristiky
figure
subplot(2,1,1)%pro H1
title('Staticke charakteristiky')
plot(Uhodnoty,statNelinearH1)
hold on
plot(Uhodnoty,statLin_1H1-0.6)
hold on
plot(Uhodnoty,statLin_2H1-0.61)
title('Staticke charakteristiky pro H1')
xlabel('vstup [m^3s^{-1}]')
ylabel('H1 [m]')
grid on
legend('nelinearni model','linearizovany model', 'linearizovany mode pro 1.2*Q1','Location','northwest')

subplot(2,1,2)%pro H2
plot(Uhodnoty,statNelinearH2)
hold on
plot(Uhodnoty,statLin_1H2-0.4)
hold on
plot(Uhodnoty,statLin_2H2-0.4)
title('Staticke charakteristiky pro H2')
xlabel('vstup [m^3s^{-1}]')
ylabel('H2 [m]')
grid on
legend('nelinearni model','linearizovany model', 'linearizovany mode pro 1.2*Q1','Location','northwest')

%% vypocteni prenosu ze stavoveho popisu
p=tf('p');
A=A_1;
B=B_1;
C=C_1;
D=D_1;
%pouziti vztahu pro vypocet prenosu
Fp=C*inv(p*eye(2)-A)*B;

%% impulzni funkce pro H2 - charakteristika, porovnani
%impulzni charakteristika ziskana z linearizovaneho modelu
sysLin2=ss(A,B,[0,1],0);%system pouze pro H2, proto matice C je [0,1]
figure
impulse(sysLin2)
hold on
syms t
impF=194.0285*(-exp(-0.3421*t)+exp(-0.0329*t));
fplot(impF,[0 250],'r')
legend('linearizovany model','impulsni funkce')
title('Impulsni charakteristika')
xlabel('cas [s]')
ylabel('H2 [m]')

%% prechodova funkce pro H2 - charakteristika, porovnani
sysLin2=ss(A,B,[0,1],0);
figure
step(sysLin2);
hold on
syms t;
prechF=-(194.0285/0.3421)*(1-exp(-0.3421*t))+(194.0285/0.0329)*(1-exp(-0.0329*t));
fplot(prechF,[0 250],'r')
legend('linearizovany model','prechodova funkce','Location','southeast')
title('Prechodova charakteristika')
xlabel('cas [s]')
ylabel('H2 [m]')

%% e

%% stavovy model frobeniova reprezentace - vypocteno na papire
Af=[0,1;-0.01125,-0.375];
Bf=[0;1];
Cf=[60,0];
Df=0;
sysLinFrob=ss(Af,Bf,Cf,Df);

%% stavovy model Jordanuv tvar
%najdu vlastni cisla matice A naseho puvodniho linearizovaneho systemu
[vlVektory,Aspektralni]=eig(A_1);
%sestavim z nich spektralni matici coz  je i nase nova matice Aj
Aj=Aspektralni;
%zapisi transformacni matici Tj a jeji inverzi
Tj=inv(vlVektory);
Tjinv=vlVektory;
%dopoctu matice Bj, Cj a Dj
Bj=Tj*B;
Cj=[0,1]*Tjinv;
Dj=0; %matice D=0 puvodniho systemu je rovna nove matici Dj
sysLinJordan=ss(Aj,Bj,Cj,Dj);

%% porovnani stavovych reprezentaci
%prechodove charakteristiky
figure
step(sysLin2,250)%prenos puvodni stavove reprezentace
hold on
step(sysLinFrob,250)%Frobeniova stavova reprezentace systemu
hold on
step(sysLinJordan,250)%Jordanova stavova reprezentace systemu
title('Prechodove charakteristiky')
xlabel('cas [s]')
ylabel('H2 [m]')
legend('puvodni system','Frobeniuv tvar','Jordanuv tvar','Location','southeast')

%impulsni charakteristiky
figure
impulse(sysLin2)%prenos puvodni stavove reprezentace
hold on
impulse(sysLinFrob)%Frobeniova stavova reprezentace systemu
hold on
impulse(sysLinJordan)%Jordanova stavova reprezentace systemu
title('Impulsni charakteristiky')
xlabel('cas [s]')
ylabel('H2 [m]')
legend('puvodni system','Frobeniuv tvar','Jordanuv tvar')

%% riditelnost a pozorovatelnost
%puvodni reprezentace
rank(A_1)
Qr_p=[B_1,A_1*B_1];
rank(Qr_p)
Qp_p=[0,1;[0,1]*A_1];
rank(Qp_p)
%Frobeniova reprezentace
rank(Af)
Qr_f=[Bf,Af*Bf];
rank(Qr_f)
Qp_f=[Cf;Cf*Af];
rank(Qp_f)
%Jordanova reprezentace
rank(Aj)
Qr_j=[Bj,Aj*Bj];
rank(Qr_j)
Qp_j=[Cj;Cj*Aj];
rank(Qp_j)
%% druha cast %%

%% a) vypocteni noveho stavoveho systemu na papire s nasledujicimi vysledky
prenosCerpadla=(3*10^(-4))/(p+2);
A_3=[-0.15,0.15,400;0.15,-0.225,0;0,0,-2];
B_3=[0;0;3*10^(-4)];
C_3=[1,0,0;0,1,0];
D_3=[0;0];

sys_v2_0=ss(A_3,B_3,C_3,D_3);%system s cerpadlem

prenosy=C_3*inv(p*eye(3)-A_3)*B_3;%prenosove funkce pro hladiny H1 a H2

%% b)
%experimentalni urceni parametru K a Ti
Fs=prenosy(2);%prenos od napeti na kotve motoru na vysku hladiny H2
K=4;
Ti=50;
Fr=K*(1+(1/(Ti*p)));%prenos PI regulatoru
Fo=Fs*Fr;
Fu=(Fo)/(1+Fo);%prenos uzavreneho systemu
%step(Fu)
%figure
%impulse(Fu)
roots([50 118.8 38.06 1.125])
%% porovnani neregulovane a regulovane soustavy

figure
step(Fs)%neregulovana soustava
hold on
step(Fu)%regulovana soustava
ylabel('H2 [m]')
xlabel('cas')
title('Porovnani prechodovych charakteristik')
legend('neregulovana soustava','regulovana soustava','Location','southeast')

dveNadrzeRizeni=sim('dve_nadrze_rizeni');
figure
plot(dveNadrzeRizeni.simout.Time,dveNadrzeRizeni.simout.Data(:,1))
title('Vystup rizeni systemu')
ylabel('Uk [V]')
xlabel('cas [s]')
%% body frekvencni charakteristiky
freq=0.5;
time=100;
frekvChar=sim('dve_nadrze_rizeni_frekvencni');
figure
plot(frekvChar.simout.Time,frekvChar.simout.Data(:,1))
hold on
plot(frekvChar.simout.Time,frekvChar.simout.Data(:,2))
xlabel('cas [s]')
ylabel('vystup [m]')
title('Frekvenci odezva na frekvenci 0.5 rad/s')
legend('vystup ot. systemu','harm. vstup','Location','northwest')

%% srovnani Nyquistovy krivky s vypoctenymi body
nyquist(Fo)
hold on
plot(0.892278,-6.783014,'Marker','+','color','red','lineWidth',1.5,'markerSize',10)
hold on
plot(0.029588,-1.872966,'Marker','+','color','yellow','lineWidth',1.5,'markerSize',10)
hold on
plot(-0.231169,-0.947920,'Marker','+','color','magenta','lineWidth',1.5,'markerSize',10)
hold on
plot(-0.108843,-0.035658,'Marker','+','color','cyan','lineWidth',1.5,'markerSize',10)
legend('Nyquistova krivka','Bod 1','Bod 2','Bod 3','Bod 4','Location','northwest')

%% d)

Kkrit=(19*2793494^(1/2))/120 - 1322/5;
%0.234356497581;
Ti=0.5;
Fr=Kkrit*(1+(1/(Ti*p)));%prenos PI regulatoru
Fo=Fs*Fr;
Fu=Fo/(1+Fo);
step(Fu,10000)
ylabel('H2 [m]')

%% e)
clc
Ti=1/0.03745;
K=3.7525;
Tdz=0;%12.518;%dopravni zpozdeni
Fdz=exp(-Tdz*p);
Fr=K*(1+(1/(Ti*p)));
Fs=prenosy(2);%prenos od napeti na kotve motoru na vysku hladiny H2
Fco=Fdz*Fs*Fr;
Fcu=Fco/(1+Fco);
step(Fcu)
hold on

Ti=50;%1/0.03745;%50;
K=4;%3.7525;%4;
Fr=K*(1+(1/(Ti*p)));
Fo=Fs*Fr;
Fu=1/(1/Fo+1);%Fo/(1+Fo);
step(Fu)
title('Odezva na jednotkový skok')
ylabel('H2 [m]')
legend('Regulátor naladěný pomocí GMK','Regulátor naladěný experimentálně','Location','southeast')

%figure
%bode(Fcu)
