%%Lipid Data 


close all
clear all
%  PWrs2 = 0;%PWrs/(10);%PWrs2/100;
  

depth = 8160 % sample sigmnals
Temp1 = 4; % cooled temperature
Temp2 = 19;% high temperature

% [dosmooth dosmooth_H2O] = [0 0],[0 2], [1 1],[1 3]
dosmooth = 1% I chnaged the formular for 1 ;0 = normalize 1= normalize +smooth + 2 = V2(unnomalized water 4C) normalized 3 = smooth+2
dosmooth_H2O = 1
Min_Ht2ORT = [];
Max_Ht2ORT = [];


      
 load('G:\Experiments\Archive2018\Temp\Lipid\Lipid_info')
    actualPowersInit = mean(PWrs2');    
    [Mat_p2p_Lipid,Mat_hlp_Lipid,SD_p2p_Lipid,SNR_p2p_Lipid] = ...
        SD_continuous_ch1(spectra,SelectedWLs,PWrs2,actualPowersInit,...
    1250,1600,1,0,depth,0); % analysis Function
    Mat_p2p_Lipid(:,2:end) = Mat_p2p_Lipid(:,2:end); % delete Wavelength info from data
    Trusted = 1:spectraCount;% selectiv plot : do not chnage(all are selected)
    L = length(SelectedWLs)*spectraCount;
    Mat_p2p= Mat_p2p_Lipid; % Just rename of peak to peak matrix
    
    a = 9; %% WL = 930 nm
    b= 17; %WL = 970=0 nm
    Division = double((Mat_p2p(a,Trusted+1)./Mat_p2p(b,Trusted+1))); % division of 930nm/970nm

        Time_Rounds=(reshape(timePoint(1:L),length(SelectedWLs),spectraCount)); %% rows+ WL; Columns= rounds
        SelectedTemp = [2.5 4 5.2 7.8 9.9 12.1 13.8 16 17.8 19]; % Temperatures to be plotted
 

        %% plot
        
        Ref_sec = (MinSecVal(:,1)*60+MinSecVal(:,2));
 [~, ind] = unique(Ref_sec)
 duplicate_ind = setdiff(1:size(Ref_sec, 1), ind);
 Ref_sec(duplicate_ind) = [];
 MinSecVal(duplicate_ind,:)=[];
 
 
Interp_TempVal = interp1(Ref_sec,MinSecVal(:,3),Time_Rounds(1,:),'spline'); 

[val,loc1] = min(abs(Interp_TempVal-4));
[val,loc2] = min(abs(Interp_TempVal-18.9));

%  
% plot(Ref_sec,MinSecVal(:,3))
% hold on


    Spectra_4C = smooth((Mat_p2p(:,loc1+1)-min(Mat_p2p(:,loc1+1)))/(max(Mat_p2p(:,loc1+1))-min(Mat_p2p(:,loc1+1))));
    Spectra_end = smooth((Mat_p2p(:,loc2+1)-min(Mat_p2p(:,loc2)))/(max(Mat_p2p(:,loc2+1))-min(Mat_p2p(:,loc2+1))));

figure(4)
subplot(2,3,2)
hold on;plot(Mat_p2p(:,1),Spectra_4C)
hold on;plot(Mat_p2p(:,1),Spectra_end)

wl_lipid = Mat_p2p(:,1)
Spectra_4C_lipid = Spectra_4C ;
Spectra_end_lipid = Spectra_end;
legend('lipid  4 C', 'lipid 19 C');
        %%
        
% add text (TEMP) to plot
l = 1;
% SlectedDegree = [1 14 52 62 70 78 87 98 108 118 127 137 147 157 167 177 187 197 207]; % location of interested temperatures in time matrix
for c = 1:length(SelectedTemp)
%    i = SlectedDegree(c)
   [valTemp,locTemp] = min(abs(MinSecVal(:,3)- SelectedTemp(c)));
    Dist = sqrt(Time_Rounds(b,:).^2 - Timeinsec(locTemp)^2);
[val, loc] = min(Dist);
locs(l) = loc; % where is the minimum Dis
locx(l) = (Time_Rounds(b,loc));
locy(l) = (Division(loc));
Temp(l) = MinSecVal(locTemp,3)

hold on
figure(4)
 subplot(2,3,1)
% to put text on plot of time-Grad
Hlb = double(Mat_p2p(:,loc+1));
 plot(Mat_p2p(:,1),Hlb)
 str = {'\leftarrow T = ',num2str(MinSecVal(locTemp,3)),'C'}
 L=double(Mat_p2p(end-5,1));
%  newStr = join(str,1);
%  text(L,Hlb(end-5),str)
l =l+1;
xlabel('Wavelength (nm)')
ylabel('I_O_A (a.u)')
end

 figure(4) % single WL or division of wls 
 subplot(2,3,3)
 g= plot(SelectedTemp,(Division(locs)),'-b');
 hold on
 set( g,'LineWidth',2);
 set(gca,'FontSize',12);
 xlim([SelectedTemp(1) SelectedTemp(end)])
 xlabel('Time(sec)');
 ylabel('Grad 1100(nm)/1000(nm)');


%   figure
  
