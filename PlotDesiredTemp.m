function [Min_Ht2ORT,Max_Ht2ORT] = PlotDesiredTemp(MinSecVal,Time_Rounds,Mat_p2p,Temp1,Temp2,dosmooth,Min_Ht2ORT,Max_Ht2ORT)
     % MinSecVal ; 1st and second columns are time (minutes & second), the third value is thermometer valeus
     %Time_Rounds ; is time we recorded OA measurments at each wavelength where 23(number of wavelength) X 74 (we repeated measurment 74 times)
     % Mat_p2p: OA peak to peak values
     % Temp1 = first interested wavelength which we want to plot
     % Temp2 = second interested wavelength which we want to plot
     
     
                         Ref_sec = (MinSecVal(:,1)*60+MinSecVal(:,2));
                         [~, ind] = unique(Ref_sec)
                         duplicate_ind = setdiff(1:size(Ref_sec, 1), ind);
                         Ref_sec(duplicate_ind) = [];
                         MinSecVal(duplicate_ind,:)=[];

                        Interp_TempVal = interp1(Ref_sec,MinSecVal(:,3),Time_Rounds(1,:),'spline'); 

                        [val,loc1] = min(abs(Interp_TempVal-Temp1));
                        [val,loc2] = min(abs(Interp_TempVal-Temp2));

                        %  
                        % plot(Ref_sec,MinSecVal(:,3))
                        % hold on
                         
                        switch dosmooth
                            case 0 % normalize
                            Spectra_4C = ((Mat_p2p(:,loc1+1)-min(Mat_p2p(:,loc1+1)))/(max(Mat_p2p(:,loc1+1))-min(Mat_p2p(:,loc1+1))));
                            Spectra_19C = (Mat_p2p(:,loc2+1)-min(Mat_p2p(:,loc2)))/(max(Mat_p2p(:,loc2+1))-min(Mat_p2p(:,loc2+1)));
                            
                            case 1 % normalize+smooth
                            Spectra_4C = ((Mat_p2p(:,loc1+1)-min(Mat_p2p(:,loc1+1)))/(max(Mat_p2p(:,loc1+1))-min(Mat_p2p(:,loc1+1))));
                            Spectra_19C = ((Mat_p2p(:,loc2+1)-min(Mat_p2p(:,loc2)))/(max(Mat_p2p(:,loc2+1))-min(Mat_p2p(:,loc2+1))));
                            Spectra_4C_smooth = smooth((Mat_p2p(:,loc1+1)-min(Mat_p2p(:,loc1+1)))/(max(Mat_p2p(:,loc1+1))-min(Mat_p2p(:,loc1+1))));
                            Spectra_19C_smooth = smooth((Mat_p2p(:,loc2+1)-min(Mat_p2p(:,loc2)))/(max(Mat_p2p(:,loc2+1))-min(Mat_p2p(:,loc2+1))));
                            case 2 % V2 normalize
                            Spectra_4C = (Mat_p2p(:,loc1+1))%/norm(Mat_p2p(:,loc1+1)));
                            Spectra_19C = (Mat_p2p(:,loc2+1)-min(Mat_p2p(:,loc2)))/(max(Mat_p2p(:,loc2+1))-min(Mat_p2p(:,loc2+1)))
                            case 3 % V2  4C unnormalize +smooth
                            Spectra_4C = smooth(Mat_p2p(:,loc1+1))%/(max(Mat_p2p(:,loc1+1))-min(Mat_p2p(:,loc1+1))));
                            Spectra_19C = smooth(Mat_p2p(:,loc2+1));
                            case 4
                           Spectra_4C = Mat_p2p(:,loc1+1);
                            Spectra_19C = Mat_p2p(:,loc2+1);
                            
                            case 5 % if water + we want to normalize to maxim of water
                             Max_Ht2ORT = max(Mat_p2p(:,loc2+1));
                             Min_Ht2ORT = min(Mat_p2p(:,loc2+1));   
                                   
                            Spectra_19C = smooth((Mat_p2p(:,loc2+1))/(Max_Ht2ORT));
                            Spectra_4C = smooth((Mat_p2p(:,loc1+1))/(Max_Ht2ORT));
                            case 6
%                              Max_Ht2ORT = max(Mat_p2p(:,loc2+1));
%                              Min_Ht2ORT = min(Mat_p2p(:,loc2+1));   
                                
                            Spectra_19C = smooth((Mat_p2p(:,loc2+1))/(Max_Ht2ORT));
                            Spectra_4C = smooth((Mat_p2p(:,loc1+1))/(Max_Ht2ORT));
                        end

%                         hold on;plot(Mat_p2p(:,1),Spectra_4C)
%                         hold on;plot(Mat_p2p(:,1),Spectra_end,'k--')
%                         
                           set(findall(gca, 'Type', 'Line'),'LineWidth',2);
                          
% p1h=plot(Mat_p2p(:,1),Spectra_19C,'k--',Mat_p2p(:,1),Spectra_4C,'LineWidth',1)
hold on
 p1h=plot(Mat_p2p(:,1),Spectra_19C_smooth,'d-',Mat_p2p(:,1),Spectra_4C_smooth,'d-','LineWidth',1)
l1h=legend('19C','4C')

xlabel('\lambda, nm');
ylabel('I_O_A (A.U)');

                         
end
                        