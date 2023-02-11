clc
clear

dataSetRTEx=[1];
dataSetACC=[1];
dataSetRTStd=[1];
hgf2=[1];lawson=[1];
for subid=0:50


    result6Hz = load(['sub',num2str(subid),'/result6Hz.mat']).result6Hz;
    result10Hz = load(['sub',num2str(subid),'/result10Hz.mat']).result10Hz;
    resultSham = load(['sub',num2str(subid),'/resultSham.mat']).resultSham;
    results=[result6Hz,result10Hz,resultSham];
    for i=1:3             
        result=results(i);
        RT = result.response.reactionTime;
        noise = result.stimulate.noiseLevel;
        stiView=result.stimulate.isHouse;
        stiSound=result.stimulate.isHighSound;
        stiProb=result.config.probOfHouse_highSound;
        ex = ((stiProb==0.9).*((stiSound==1).*(stiView==1)+(stiSound==0).*(stiView==0))+...
        (stiProb==0.1).*((stiSound==0).*(stiView==1)+(stiSound==1).*(stiView==0)))+...
        ((stiProb==0.9).*((stiSound==1).*(stiView==0)+(stiSound==0).*(stiView==1))+...
        (stiProb==0.1).*((stiSound==1).*(stiView==1)+(stiSound==0).*(stiView==0)))*-1;
        c=double(result.response.isLeftKeyDown==result.stimulate.isHouse);

        c(c==-1)=NaN;
        RT(RT<0)=NaN;

        RTmodify=RT(:,:);
        RTmodify(c~=1)=NaN;
        RTSTD=nanstd(RTmodify);
        RTMean=nanmean(RTmodify)  ;  
        RTmodify(RTmodify>RTMean+3*RTSTD)=NaN;
        RTmodify(RTmodify<RTMean-3*RTSTD)=NaN;

        RTmodify=1000.*RTmodify;
        
        dataSetRTEx(subid+1,1)=subid;
        dataSetACC(subid+1,1)=subid;
        dataSetRTStd(subid+1,1)=subid;

        dataSetRTEx(subid+1,1+(i-1)*9+1)=nanmean(RTmodify((noise==1)+(ex==1)==2));
        dataSetRTEx(subid+1,1+(i-1)*9+2)=nanmean(RTmodify((noise==2)+(ex==1)==2));
        dataSetRTEx(subid+1,1+(i-1)*9+3)=nanmean(RTmodify((noise==3)+(ex==1)==2));
        dataSetRTEx(subid+1,1+(i-1)*9+4)=nanmean(RTmodify((noise==1)+(ex==-1)==2));
        dataSetRTEx(subid+1,1+(i-1)*9+5)=nanmean(RTmodify((noise==2)+(ex==-1)==2));
        dataSetRTEx(subid+1,1+(i-1)*9+6)=nanmean(RTmodify((noise==3)+(ex==-1)==2));
        dataSetRTEx(subid+1,1+(i-1)*9+7)=nanmean(RTmodify);
        dataSetRTEx(subid+1,1+(i-1)*9+8)=result.response.correctRace;
        dataSetRTEx(subid+1,1+(i-1)*9+9)=lillietest(RTmodify);

        dataSetACC(subid+1,1+(i-1)*7+1)=nanmean(c((noise==1)+(ex==1)==2));
        dataSetACC(subid+1,1+(i-1)*7+2)=nanmean(c((noise==2)+(ex==1)==2));
        dataSetACC(subid+1,1+(i-1)*7+3)=nanmean(c((noise==3)+(ex==1)==2));
        dataSetACC(subid+1,1+(i-1)*7+4)=nanmean(c((noise==1)+(ex==-1)==2));
        dataSetACC(subid+1,1+(i-1)*7+5)=nanmean(c((noise==2)+(ex==-1)==2));
        dataSetACC(subid+1,1+(i-1)*7+6)=nanmean(c((noise==3)+(ex==-1)==2));
        dataSetACC(subid+1,1+(i-1)*7+7)=result.response.correctRace;



        dataSetRTStd(subid+1,1+(i-1)*7+1)=nanstd(RT((noise==1)+(ex==1)==2));
        dataSetRTStd(subid+1,1+(i-1)*7+2)=nanstd(RT((noise==2)+(ex==1)==2));
        dataSetRTStd(subid+1,1+(i-1)*7+3)=nanstd(RT((noise==3)+(ex==1)==2));
        dataSetRTStd(subid+1,1+(i-1)*7+4)=nanstd(RT((noise==1)+(ex==-1)==2));
        dataSetRTStd(subid+1,1+(i-1)*7+5)=nanstd(RT((noise==2)+(ex==-1)==2));
        dataSetRTStd(subid+1,1+(i-1)*7+6)=nanstd(RT((noise==3)+(ex==-1)==2));
        dataSetRTStd(subid+1,1+(i-1)*7+7)=nanstd(RT);

        RT = result.response.reactionTime;
         RT(RT<0)=NaN;
        RT(c~=1) =NaN;
        RTSTD=nanstd(RT);
        RTMean=nanmean(RT)  ;  
        RT(RT>RTMean+3*RTSTD)=NaN;
        RT(RT<RTMean-3*RTSTD)=NaN;
        logRT=log(1000.*RT);

        input = double(result.stimulate.isHouseFace_highLowSound);
        uncertainty = result.stimulate.noiseUncertainty;

        temp2= tapas_fitModel(logRT',...
             [input;uncertainty]',...
             'my_tapas_hgf_binary_pu_tbt_config',...
             'my_tapas_logrt_linear_binary_3_config',...
             'tapas_quasinewton_optim_config');
    
        lawson(subid+1,i)=temp2.p_prc.om(2);
        lawson(subid+1,i+3)=temp2.p_prc.om(3);
    
        lawson(subid+1,i+21)=temp2.optim.LME;

        lawson(subid+1,i+6)=temp2.p_obs.be0;
        lawson(subid+1,i+9)=temp2.p_obs.be1;
        lawson(subid+1,i+12)=temp2.p_obs.be2;
        lawson(subid+1,i+15)=temp2.p_obs.be3;
        lawson(subid+1,i+18)=temp2.p_obs.ze;
        
        subinput = double(result.response.isCorrect);
        for j=1:320
            if(subinput(j)==0)
                subinput(j)=1-input(j);
            elseif(subinput(j)==1)
                subinput(j)=input(j);
            else
                subinput(j)=NaN;
            end
        end
        tempRW= tapas_fitModel(subinput',...
             input',...
             'tapas_rw_binary_config',...
             'tapas_unitsq_sgm_config',...
             'tapas_quasinewton_optim_config');
        tempSK1= tapas_fitModel(subinput',...
             [input]',...
             'tapas_sutton_k1_binary_config',...
             'tapas_unitsq_sgm_config',...
             'tapas_quasinewton_optim_config');
        tempHGF2= tapas_fitModel(logRT',...
             [input;uncertainty]',...
             'my_tapas_hgf_binary_pu_tbt_config',...
             'my_tapas_logrt_linear_binary_config',...
             'tapas_quasinewton_optim_config');


        hgf2(subid+1,i)=tempHGF2.p_prc.om(2);
        hgf2(subid+1,i+3)=tempHGF2.p_prc.om(3);
    
        hgf2(subid+1,i+21)=tempHGF2.optim.LME;

        hgf2(subid+1,i+6)=tempHGF2.p_obs.be0;
        hgf2(subid+1,i+9)=tempHGF2.p_obs.be1;
        hgf2(subid+1,i+12)=tempHGF2.p_obs.be2;
        hgf2(subid+1,i+15)=tempHGF2.p_obs.be3;
        hgf2(subid+1,i+18)=tempHGF2.p_obs.ze;

        if(i==1)
            tACS = '6Hz';
        elseif i==2
            tACS = '10Hz';
        else
            tACS = 'Sham';
        end

         dataLawson=temp2;

        save(['model/modelLawson','Sub',num2str(subid),'tACS',tACS,'.mat'],'dataLawson')
         dataRW=tempRW;
        save(['model/modelRW','Sub',num2str(subid),'tACS',tACS,'.mat'],'dataRW')
        dataSK1=tempSK1;
        save(['model/modelSK1','Sub',num2str(subid),'tACS',tACS,'.mat'],'dataSK1')
        dataHGF2=tempHGF2;
        save(['model/modelHGF2','Sub',num2str(subid),'tACS',tACS,'.mat'],'dataHGF2')

%% 下面是作图时留下的代码
%         figure()
%         subplot(3,1,1)
%         x=1:320;
%         axis([0,320,-0.3,2]);
%         hold on
%         plot(x,temp2.traj.wt(:,1),"m")
%         plot(x,result.config.probOfHouse_highSound,"black")
%         plot(x,temp2.traj.muhat(:,1),"Color",[124,194,204]/255)
%         scatter(x,input,".","blue")
    %     
    %     scatter(x,tempOutput,".","Color",[233,145,108]/255)
    %     title("logRT")
    %     legend("learningRate","P(face|highSound)","posteriorBelief","stimulus","subject'sInput")
    %     
    %     subplot(3,1,2)
    %     x=1:320;
    %     axis([0,320,-0.3,2]);
    %     hold on
    %     plot(x,temp2.traj.wt(:,1),"m")
    %     plot(x,result.config.probOfHouse_highSound,"black")
    %     plot(x,temp2.traj.muhat(:,1),"Color",[124,194,204]/255)
    %     scatter(x,input,".","blue")
    %     
    %     scatter(x,tempOutput,".","Color",[233,145,108]/255)
    %     title("unitsq")
    %     legend("learningRate","P(face|highSound)","posteriorBelief","stimulus","subject'sInput")
    % 
    % 
    %     subplot(3,1,3)
    %     hold on 
    %     plot(x,result.config.probOfHouse_highSound,"black")
    %     plot(x,result.response.reactionTime)
    %     legend('RT')
    
    end
end
%%
input = double(result6Hz.stimulate.isHouseFace_highLowSound);
subinput = double(result6Hz.response.isCorrect);
%        tempOutput = double(result.response.isCorrect);
        for j=1:320
            if(subinput(j)==0)
%                tempOutput(j)=1-input(j);
                subinput(j)=1-input(j);
            elseif(subinput(j)==1)
%                tempOutput(j)=input(j);
                subinput(j)=input(j);
            else
%                tempOutput(j)=NaN;
                subinput(j)=NaN;
            end
        end