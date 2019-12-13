%% in this code you can calculate: p2p in datail, and SD

 function [Mat_p2p,Mat_hlp,SD_p2p,SNR_p2p] = SD_continuous_ch1(spectra,selectedWLs,PWrs2,actualPowersInit,Beg,End,flagPow,Flagplot,depth,R)
%%
%  L =length( (spectra{1,1}{1}))
 Rounds = round(size(spectra,2)/(numel(selectedWLs)*depth)); % should be changed depth
%  spectra(1)=[];
spectra2=reshape(spectra,numel(spectra)/Rounds,Rounds);

for i =1:size(spectra2,2)
spectra3{i}=reshape(spectra2(:,i),numel(spectra2(:,i))/length(selectedWLs),length(selectedWLs));
end
 %%
 if R==0
 
  Rounds = size(spectra3,2);
 else
     
  Rounds = R
 end
 
  S_raw = spectra3;
 L =length( (S_raw{1,1}));
 
 

 
 
for rounds = 1:Rounds%: size(S_raw,1)
    hlbSpectr = 0;
for WL = 1:length(selectedWLs)
    
    data(WL,:)= (S_raw{rounds}(:,WL))';     
   
    hlbData = abs(hilbert(data(WL,:)));
    
    p2p(WL,:) = max(data(WL,Beg:End))-min(data(WL,Beg:End));
    singalSum = sum(hlbData(Beg:End))/(End-Beg +1); %/actualPowersInit(WL);
    hlbSpectr = [hlbSpectr singalSum ];   
   
      if Flagplot == 1  & selectedWLs(WL)== 1000
   figure(2)
    hold on; plot( data(WL,:));
  
    title(['WL = ' num2str(selectedWLs(WL))])
      else if Flagplot == 1  & selectedWLs(WL)== 930
   figure(3)
    hold on; plot( data(WL,:));
     title(['WL = ' num2str(selectedWLs(WL))])
          else if Flagplot == 1  & selectedWLs(WL)== 970
                  figure(4)
    hold on; plot( data(WL,:));
     title(['WL = ' num2str(selectedWLs(WL))])
              end
          end
      end     
    
    if Flagplot == 2   &  ~mod(WL,5) 
   figure(2)
    hold on; plot( data(WL,:));
    end
    
end
    
 p2pspectra(rounds,:) =  p2p;
 hlbSpectra(rounds,:) =  hlbSpectr;
end

hlbSpectra (:,1)='';

if flagPow == 1
       for s = 1:Rounds
              actualPowersInit = mean(PWrs2');
         hlbSpectra(s,:) = hlbSpectra(s,:) ./actualPowersInit;
          p2pspectra(s,:)  = p2pspectra(s,:)./actualPowersInit;
       end
end
       
     p2pspectra2 = single(p2pspectra);   
     hlbSpectra2 = single(hlbSpectra);
  

% selectedWLs = single(selectedWLs);


for s = 1:Rounds
    Mat = [selectedWLs;p2pspectra2(s,:)];
    Mat =  sortrows(Mat',1);
    Mat_p2p(:,1) = Mat(:,1);
    Mat_p2p(:,s+1) = Mat(:,2);
    
    clear Mat
    
     Mat = [selectedWLs;hlbSpectra2(s,:)];
     Mat = sortrows(Mat',1);
     Mat_hlp(:,1) = Mat(:,1);
     Mat_hlp(:,s+1) = Mat(:,2); 
    
  
    SD_p2p(s,:) = std(Mat(:,2:end)');
%      SNR_p2p2(round,:) = mean(Mat_p2p(:,2:end)'./max(Mat_p2p(:,2:end)'))./std(Mat_p2p(:,2:end)'./max(Mat_p2p(:,2:end)'));
   SNR_p2p(s,:) = mean(Mat(:,2:end)')./std(Mat(:,2:end)');
end  

% sig = mean(Mat_p2p(:,2:end)');
%     figure(2)
%     plot(Mat_p2p(:,1),sig)
%        
  
%     SD_hlb(round,:) = std(Mat_hlb(:,2:end)');
%     SNR_hlb(round,:) = mean(Mat_hlb(:,2:end)')./std(Mat_hlb(:,2:end)');
%     
%     sig = mean(Mat_hlb(:,2:end)');
    
end




