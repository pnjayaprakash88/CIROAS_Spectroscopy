% load and prepare
clc;
clear all;
tic;
close all;
sqAll = 0;
concentration = [0	53.33333333	103.1111111	149.5703704	306.265679	452.5146337];
corrCoef = [16.17679699; 9.021466244; 9.083708805; 21.3906428; 11.86521891; 12.67411198; 12.6326544; 12.64379383; 13.83558856; 13.9001687; 13.50744128; 14.62109247; 15.48062414; 16.88663723; 14.30898727; 15.43750824; 10.01368437; 19.50056876; 17.57366669; 31.23298398; 46.97980891; 50.91666588; 44.11586672; 40.58415353; 51.4706636; 49.58738222; 42.81094106; 30.63597302; 31.14372988; 27.04816778; 27.06289606; 22.92118457; 19.82839016; 18.92304619; 19.25264232; 18.31001255; 17.14796875; 20.19463554; 18.54021376; 17.28782905; 16.4161642; 18.97515616; 26.24242614; 25.89885351; 24.08720268; 152.1607402];
x1= concentration';
X = [ones(length(x1),1) x1];
Fs = 5e5;
flt1 = 1e2;
flt2 = 5e4;
S = 0;
flags = [0 1 2 3 5 10];
Inters = 0;
counts = 0;
bkgBeg = 750;
bkgEnd = 1300;
omits = 1:2;
omit2s =0:1;

allCounts = counts; counts=0;
for iiLoad = 1:6
    flag = flags(iiLoad);
    switch flag
        case 0
            load('20170625WaterOnly.mat');
        case 1
            load('20170625Glu1.mat');
        case 2
            load('20170625Glu2.mat');
        case 3
            load('20170625Glu3.mat');
        case 5
            load('20170625Glu5.mat');
        case 10
            load('20170625Glu10.mat');
    end
    nWL = length(SelectedWLs);
    nT = spectraCount;
    signalN = length(spectra)/nT/nWL;
    spectra = reshape ( spectra,signalN,nWL,nT);
    spectraProcessed = spectra;
    SelectedWLsSorted = sort  (SelectedWLs)';
    for ii = 1:nWL
        sortedWL = ii;
        for jj = 1:nWL
            if SelectedWLs(ii)== SelectedWLsSorted(jj)
                sortedWL = jj;
            end
        end
        spectraProcessed (:,sortedWL,:)=spectra(:,ii,:);
    end
    for ii = 1:nWL
        for jj = 1:nT
            spectraProcessed (:,ii,jj) = filtData(spectraProcessed (:,ii,jj),Fs,[flt1 flt2]);
            spectraProcessed (:,ii,jj) = abs(hilbert(spectraProcessed (:,ii,jj)));
        end
    end
    if iiLoad ==1
        spectraArray = {spectraProcessed};
    else
        spectraArray = [spectraArray {spectraProcessed}];
    end
    if iiLoad ==1
        tAllArray = {tAll};
    else
        tAllArray = [tAllArray {tAll}];
    end
end

intercepts = zeros(nWL,1);
clear 'timePoint' 'times';
clear 'sortedWL' 'powers' 'resValues';
toc

%% Parameter Optimisation
clc; close all;
count = 0;

%[m,n]=size (sSelect);
for kk = 1:2
    
    omits = 3;
    omit2s =0;
    
    b11 = 1515;
    b12 = 1600;
    e11 = 2370; % (b1 + )
    e12 = 2400;
    b22 = 1590;
    e21 = 2815;
    e22 = 2835;
    step = 5;

    
    if kk==1
        beG1 = b11;
        beG2 = b22;
        enD1 = e11;
        enD2 = e22;
    else
        beG1 = b12;
        beG2 = b22+20;
        enD1 = e12;
        enD2 = e21;
    end
    
    omit = omits;
    omit2 = omit2s;

    for iiFlag = 1:6
        flag = flags(iiFlag);
        spectraProcessed = cell2mat(spectraArray(iiFlag));
        [signalN,nWL,nT]=size (spectraProcessed);
        tAll=cell2mat(tAllArray(iiFlag)) -0.39;
        tMin = min (tAll);
        tMax = max (tAll);
        results = squeeze ( spectraProcessed(1,:,:));
        for jj = 1:nT
            for ii = 1:nWL
                x = (tAll(jj) - tMin)/(tMax-tMin);
                beGNow = round(beG2- (beG2 -beG1)*x);
                enDNow = round(enD2- (enD2 -enD1)*x);
                results(ii,jj) = sum (spectraProcessed(beGNow:enDNow,ii,jj))/(enDNow-beGNow+1) -sum (spectraProcessed (bkgBeg:bkgEnd,ii,jj))/(bkgEnd-bkgBeg+1) ;
                results(ii,jj) = results(ii,jj)*corrCoef(ii);
            end
        end
        results = results';

        interceptsSecondary = zeros(nWL,1);
        for ii = 1:nWL
            [minX,pos] = min (results (:,ii));
            if pos>omit&&pos<nT
                x1 = tAll(omit:pos-omit2);
                y1 = results (omit:pos-omit2,ii);
            end
            Xd = [ones(length(x1),1) x1];
            k1 = Xd\y1;
            intercepts(ii) = -k1(1)/k1(2);
        end
        if iiFlag == 1
            Inters =  intercepts;
        else
            Inters =  [Inters, intercepts];
        end
    end
    
    for ii = 20:24
        if max(Inters(ii,:))<5&&min(Inters(ii,:))>0
            Y = Inters(ii,:)';
            k = X\Y;
            
            squares = 0;
            for jj = 1:6
                squares = squares +(k(2)*concentration (jj) +k(1)-Y(jj))^2;
            end
            
            if k(2)< -0.002 && squares<0.0015
                %%% Plotting the intercept graph
                
                figure
                count=count +1;
                count
                squares
                SelectedWLsSorted(ii)
                kk
                
                plot (concentration, smooth(Inters(ii,:)'), '*-')
                %ii
            end
        end
    end
    
end
%% Visalisation of glucose spectra
beG1 = 1565; enD1 = 2375; beG2 = 1585; enD2 = 2820;

flag = 4;
spectraProcessed_Glucose = cell2mat(spectraArray(flag));
[signalN,nWL,nT]=size (spectraProcessed_Glucose);
tAll=cell2mat(tAllArray(flag));
tMin = min (tAll);
tMax = max (tAll);
results = squeeze ( spectraProcessed_Glucose(1,:,:));
for jj = 1:nT
    for ii = 1:nWL
        x = (tAll(jj) - tMin)/(tMax-tMin);
        beGNow = round(beG2- (beG2 -beG1)*x);
        enDNow = round(enD2- (enD2 -enD1)*x);
        results(ii,jj) = sum (spectraProcessed_Glucose(beGNow:enDNow,ii,jj))/(enDNow-beGNow+1) -sum (spectraProcessed (bkgBeg:bkgEnd,ii,jj))/(bkgEnd-bkgBeg+1) ;
        results(ii,jj) = results(ii,jj)*corrCoef(ii);
    end
end

corr2 = [0.787107334 1.074057623 1.140648137 0.504070786 0.93504937 1.160114887 1.141803784 1.216257379 1.165074165 1.141233747 1.219509758 1.118618413 1.180629159 1.076365391 1.228352508 1.121989225 0.997323264 1.090365711 0.97356989 0.843346739 0.770924943 0.752166822 0.839487187 0.920182101 0.886417881 0.748792725 0.820248312 0.895814234 1.037496173 1.116744933 1.024438116 1.067262043 1.124926277 1.102717051 1.06423691 1.138289283 1.207496554 1.228420745 1.113560523 1.238143298 1.147253575 1.016193317 1.031231805 0.994491163 0.970882747 0.709148857 ];
waterSpec = [0.06707 0.298 0.416 0.1953 0.295 0.6778 1.135 1.236 1.25 1.24 1.22 1.15 1.085 1.122 1.359 1.806 2.88 4.618 7.387 16.29 26.93 31.05 31.55 30.9 29.76 28.61 25.72 22.63 19.76 14.84 11.91 9.97 8.54 7.435 6.603 6.055 5.755 5.65 5.75 6.15 6.961 7.896 9.255 9.52 10.72 48.16];
correctionIndex = 2;
corr2 = results(:,correctionIndex)./waterSpec';
indexSelection = [2:2:14 17 30];
figure;
hold on;

for indexjj = 1:length(indexSelection)
    jj = indexSelection(indexjj);
    plot (SelectedWLsSorted(1:end-1), results(1:end-1,jj));
end

%% Visalisation of Water Spectra
beG1 = 1565; enD1 = 2375; beG2 = 1585; enD2 = 2820;

flag = 1;
spectraProcessed_Water = cell2mat(spectraArray(flag));
[signalN,nWL,nT]=size (spectraProcessed_Water);
tAll=cell2mat(tAllArray(flag));
tMin = min (tAll);
tMax = max (tAll);
results_Water = squeeze ( spectraProcessed_Water(1,:,:));
for jj = 1:nT
    for ii = 1:nWL
        x = (tAll(jj) - tMin)/(tMax-tMin);
        beGNow = round(beG2- (beG2 -beG1)*x);
        enDNow = round(enD2- (enD2 -enD1)*x);
        results_Water(ii,jj) = sum (spectraProcessed_Water(beGNow:enDNow,ii,jj))/(enDNow-beGNow+1) -sum (spectraProcessed (bkgBeg:bkgEnd,ii,jj))/(bkgEnd-bkgBeg+1) ;
        results_Water(ii,jj) = results_Water(ii,jj)*corrCoef(ii);
    end
end
waterSpec = [0.06707 0.298 0.416 0.1953 0.295 0.6778 1.135 1.236 1.25 1.24 1.22 1.15 1.085 1.122 1.359 1.806 2.88 4.618 7.387 16.29 26.93 31.05 31.55 30.9 29.76 28.61 25.72 22.63 19.76 14.84 11.91 9.97 8.54 7.435 6.603 6.055 5.755 5.65 5.75 6.15 6.961 7.896 9.255 9.52 10.72 48.16];
correctionIndex = 2;
corr2 = results_Water(:,correctionIndex)./waterSpec';
indexSelection = [2:2:14 17 30];
figure;
hold on;

for indexjj = 1:length(indexSelection)
    jj = indexSelection(indexjj);
    plot (SelectedWLsSorted(1:end-1), results_Water(1:end-1,jj));
end


%% Visalisation of Difference Spectra

figure; hold on;
plot(SelectedWLsSorted(1:2:end-3), smooth(abs(results(1:2:end-3,4)/max(results(1:2:end-3,4))-results_Water(1:2:end-3,26)/max(results_Water(1:2:end-3,26)))),'r*-');
plot(SelectedWLsSorted(1:2:end-3), smooth(abs((results(1:2:end-3,26))/max(results(1:2:end-3,26))-results_Water(1:2:end-3,26)/max(results_Water(1:2:end-3,26)))),'ks-');

figure; hold on;
plot(SelectedWLsSorted(1:2:end-1), (results_Water(1:2:end-1,5)/10)/norm(results_Water(1:2:end-1,5)/10),'k-');
plot(SelectedWLsSorted(1:2:end-1), results_Water(1:2:end-1,17)/10,'r-');
plot(SelectedWLsSorted(1:2:end-1), (results(1:2:end-1,14)/10)/(norm(results(1:2:end-1,14)/10)),'b-');
