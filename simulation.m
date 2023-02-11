%%
clc
clear

lawson=zeros(100,51,21);
for times=0:100
    for subid=0:50
        model6Hz = load(['./model/model','Lawson','Sub',num2str(subid),'tACS6Hz.mat']).(['data','Lawson']);
        model10Hz = load(['./model/model','Lawson','Sub',num2str(subid),'tACS10Hz.mat']).(['data','Lawson']);
        modelSham = load(['./model/model','Lawson','Sub',num2str(subid),'tACSSham.mat']).(['data','Lawson']);
        models=[model6Hz,model10Hz,modelSham];
    
        for i=1:3             
            model=models(i);
            
    
            simulationArray= tapas_simModel(model.u,...
                 'my_tapas_hgf_binary_pu_tbt',...
                 model.p_prc.p,...
                 'my_tapas_logrt_linear_binary_3',...
                 model.p_obs.p);
            
            simulationParams = tapas_fitModel(simulationArray.y, ...
                simulationArray.u, ...
                'my_tapas_hgf_binary_pu_tbt_config', ...
                'my_tapas_logrt_linear_binary_3_config', ...
                'tapas_quasinewton_optim_config');
            lawson(times,subid+1,i)=simulationParams.p_prc.om(2);
            lawson(times,subid+1,i+3)=simulationParams.p_prc.om(3);
        
    
            lawson(times,subid+1,i+6)=simulationParams.p_obs.be0;
            lawson(times,subid+1,i+9)=simulationParams.p_obs.be1;
            lawson(times,subid+1,i+12)=simulationParams.p_obs.be2;
            lawson(times,subid+1,i+15)=simulationParams.p_obs.be3;
            lawson(times,subid+1,i+18)=simulationParams.p_obs.ze;
            
        
        end
    end
end

simulationModel = struct;
simulationModel.simulation = lawson;
paramsExpetation = nanmean(lawson,1);
simulationModel.paramsExpetation =reshape(paramsExpetation,[51,21]);
%%

clc
dataSetRTEx=zeros(100,51,21);
dataSetRTStd=zeros(100,51,21);
for times=1:100
    for subid=0:50
        result6Hz = load(['sub',num2str(subid),'/result6Hz.mat']).result6Hz;
        result10Hz = load(['sub',num2str(subid),'/result10Hz.mat']).result10Hz;
        resultSham = load(['sub',num2str(subid),'/resultSham.mat']).resultSham;
        results=[result6Hz,result10Hz,resultSham];
        model6Hz = load(['./model/model','Lawson','Sub',num2str(subid),'tACS6Hz.mat']).(['data','Lawson']);
        model10Hz = load(['./model/model','Lawson','Sub',num2str(subid),'tACS10Hz.mat']).(['data','Lawson']);
        modelSham = load(['./model/model','Lawson','Sub',num2str(subid),'tACSSham.mat']).(['data','Lawson']);
        models=[model6Hz,model10Hz,modelSham];
        for i=1:3             
            result=results(i);
            noise = result.stimulate.noiseLevel;
            stiView=result.stimulate.isHouse;
            stiSound=result.stimulate.isHighSound;
            stiProb=result.config.probOfHouse_highSound;
            ex = ((stiProb==0.9).*((stiSound==1).*(stiView==1)+(stiSound==0).*(stiView==0))+...
            (stiProb==0.1).*((stiSound==0).*(stiView==1)+(stiSound==1).*(stiView==0)))+...
            ((stiProb==0.9).*((stiSound==1).*(stiView==0)+(stiSound==0).*(stiView==1))+...
            (stiProb==0.1).*((stiSound==1).*(stiView==1)+(stiSound==0).*(stiView==0)))*-1;
            c=double(result.response.isLeftKeyDown==result.stimulate.isHouse);

            model=models(i);
            
    
            simulationArray= tapas_simModel(model.u,...
                 'my_tapas_hgf_binary_pu_tbt',...
                 model.p_prc.p,...
                 'my_tapas_logrt_linear_binary_3',...
                 model.p_obs.p);
            
            simulationY = exp(simulationArray.y);
    
            dataSetRTEx(times,subid+1,(i-1)*7+1)=nanmean(simulationY((noise==1)+(ex==1)==2));
            dataSetRTEx(times,subid+1,(i-1)*7+2)=nanmean(simulationY((noise==2)+(ex==1)==2));
            dataSetRTEx(times,subid+1,(i-1)*7+3)=nanmean(simulationY((noise==3)+(ex==1)==2));
            dataSetRTEx(times,subid+1,(i-1)*7+4)=nanmean(simulationY((noise==1)+(ex==-1)==2));
            dataSetRTEx(times,subid+1,(i-1)*7+5)=nanmean(simulationY((noise==2)+(ex==-1)==2));
            dataSetRTEx(times,subid+1,(i-1)*7+6)=nanmean(simulationY((noise==3)+(ex==-1)==2));
            dataSetRTEx(times,subid+1,(i-1)*7+7)=nanmean(simulationY);
    

    
            dataSetRTStd(times,subid+1,(i-1)*7+1)=nanstd(simulationY((noise==1)+(ex==1)==2));
            dataSetRTStd(times,subid+1,(i-1)*7+2)=nanstd(simulationY((noise==2)+(ex==1)==2));
            dataSetRTStd(times,subid+1,(i-1)*7+3)=nanstd(simulationY((noise==3)+(ex==1)==2));
            dataSetRTStd(times,subid+1,(i-1)*7+4)=nanstd(simulationY((noise==1)+(ex==-1)==2));
            dataSetRTStd(times,subid+1,(i-1)*7+5)=nanstd(simulationY((noise==2)+(ex==-1)==2));
            dataSetRTStd(times,subid+1,(i-1)*7+6)=nanstd(simulationY((noise==3)+(ex==-1)==2));
            dataSetRTStd(times,subid+1,(i-1)*7+7)=nanstd(simulationY);
    
        
        
    

        
        end
    end
end
REExExpectation=reshape(nanmean(dataSetRTEx,1),[51,21]);
REStdExpectation=reshape(nanmean(dataSetRTStd,1),[51,21]);
behaviorRefcovery = struct;
behaviorRefcovery.simulation.RTEx=dataSetRTEx;
behaviorRefcovery.simulation.RTStd=dataSetRTStd;
behaviorRefcovery.behaviorExpectation.RTEx=REExExpectation;
behaviorRefcovery.behaviorExpectation.RTStd=REStdExpectation;