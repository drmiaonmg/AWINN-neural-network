clc
 
clear all
 
close all

tic
% parallel computing
delete(gcp)
parpool;
filename='data.xlsx';
sheet=1;
xlRange = 'A2:AC37';
Column1=xlsread(filename,sheet,xlRange); %查找输入
filename='data.xlsx';
sheet=1;
xlRange = 'AD2:AF37';
Column2=xlsread(filename,sheet,xlRange); %查找输入

p=Column1';
  
t=Column2(:,1)';
 
[p1,ps]=mapminmax(p,0,1);
 
[t1,ts]=mapminmax(t,0,1);

%[trainInd,valInd,testInd] = dividerand(Q,trainRatio,valRatio,testRatio)

% trainsample.p=p1(:,1:24);
% valsample.p=p1(:,27:31);
testsamplep=p1(:,31:36);
% trainsample.t=t1(:,1:26);
%valsample.t=t1(:,31:36);
testsamplet=t1(:,31:36);
% func={'trainlm','traingdx','traingda','traingd','trainbr','trainscg'};
% func1=size(func);
% trans={'tansig','logsig'};
% trans1=size(trans);


x = p1;
p2=p1;
t = t1;
t2=t1;
zuheshu=0;
% Combine data from different group to train neural network model
n=1:30;
i001=1;
for ju1=1:11
    zuhe1=[ju1 30];
    for ju2=12:13
        zuhe2=[zuhe1 ju2];
        for ju3=15:16
            zuhe3=[zuhe2 ju3];
            for ju4=17:25
                zuhe4=[zuhe3 ju4];
                for ju5=26:29
                    zuhe5=[zuhe4 ju5];
                    zuhe(i001,1:6)=zuhe5;
                    i001=i001+1;
                end
            end
        end
    end
end

for j=1:(i001-1)
    trainxulie(j,1:24)=setdiff(n,zuhe(j,:));
    trainvaluep(((j-1)*29+1):j*29,1:24)=p1(:,trainxulie(j,1:24));
     valvaluep(((j-1)*29+1):j*29,1:6)=p1(:,zuhe(j,:));
    trainvalue(j,1:24)=t1(1,trainxulie(j,1:24));
    validatevalue(j,1:6)=t1(1,zuhe(j,:));
end
% parallel computing begin
parfor j=1:(i001-1)
%     trainSize1=zeros(1,1584);
%     trainSize2=zeros(1,1584);
    pbias1=zeros(1584,24);
    pbias2=zeros(1584,6);
    pbias3=zeros(1584,6);
    erfang=zeros(1584,12);
    trainxulie1=trainxulie(j,:);
    valxulie1=zuhe(j,:);
    for i=1:20000
        trainFcn = 'trainlm';  % Bayesian Regularization
        hiddenLayerSize =20;
        hiddenLayerSize1=13;
        hiddenLayerSize2=6;
        hiddenLayerSize3=10;
        %setting the structure of neural network
        net = feedforwardnet ([hiddenLayerSize,hiddenLayerSize1,hiddenLayerSize2],trainFcn);
        %net = feedforwardnet ([hiddenLayerSize,hiddenLayerSize1,hiddenLayerSize2,hiddenLayerSize3],trainFcn);
        % Setup Division of Data for Training, Validation, Testing
        %RandStream.setGlobalStream(RandStream('mt19937ar','seed',1)); % to get constant result
        net.divideFcn = 'divideind'; % Divide targets into three sets using blocks of indices
        [trainInd,valInd,testInd] = divideind(36,trainxulie1,valxulie1,31:36);
        % net.divideParam.trainInd = 1;
        % net.divideParam.valRatio = 0;
        % net.divideParam.testRatio = 0;
        net.trainParam.showwindow=0;
        net.divideParam.trainInd=trainInd;
        net.divideParam.valInd=valInd;
        net.divideParam.testInd=testInd;
        %TRAINING PARAMETERS
        net.trainParam.show=50;  % of ephocs in display
        net.trainParam.lr=0.05;  % learning rate
        net.trainParam.max_fail=800;
        net.trainParam.epochs=1000000;  % max epochs
        net.trainParam.goal=0.00001;  % training goal
        net.performFcn='mse';  % Name of a network performance function %type help nnperformance
        %net.numLayers = 3;  
        %net.layerConnect(3,1) = 1;
        %net.outputConnect = [0 0 1 ];
        %net.layers{1}.size = 12; 
        net.layers{1}.transferFcn = 'poslin';
        %net.layers{2}.size = 20; 
        net.layers{2}.transferFcn = 'poslin';
         net.layers{3}.transferFcn = 'logsig';
%          net.layers{3}.transferFcn = 'logsig';
        %net.outputs{2}.range = [1 20];
        % Train the Network


        [net,tr] = train(net,x,t); 
        % Test the Network
        y0 = net(trainvaluep(((j-1)*29+1):j*29,:));
        e = gsubtract(trainvalue(j,1:24),y0);
        performance = perform(net,trainvalue,y0)
        y1 = net(valvaluep(((j-1)*29+1):j*29,:));
        performance1 = perform(net,validatevalue,y1)
        y2 = net(testsamplep);
        performance2 = perform(net,testsamplet,y2)
        % View the Network

        %trainoutput=mapminmax('reverse',y,ts);
        %trainvalue=mapminmax('reverse',trainsample.t,ts);
        trainoutput=y0;
        %trainvalue=trainsample.t;

        ave_obs = sum(trainvalue(j,:))/numel(trainvalue(j,:));   
        Numerator = sum(power(trainvalue(j,:)-trainoutput,2)); 
        Denominator = sum(power(trainvaluep(j,:)-ave_obs,2));
        nse1 = 1 - Numerator/Denominator;

        [r2, rmse] = rsquare(trainoutput,trainvalue(j,:));
        errors=trainvalue(j,:)-trainoutput;
        trainSize=size(trainoutput,2);
        for q=1:trainSize
            zongshu1=abs(trainvalue(j,q));
            piancha1=abs(trainvalue(j,q)-trainoutput(q));
            pbias1(j,q)=100*piancha1/zongshu1;
        end
        pbi1=max(pbias1(j,:));

        validateoutput=y1;
        %validatevalue=valsample.t;
        %validateoutput=mapminmax('reverse',y1,ts);
        %validatevalue=mapminmax('reverse',valsample.t,ts);

        ave_obs1 = sum(validatevalue(j,:))/numel(validatevalue(j,:));   
        Numerator = sum(power(validatevalue(j,:)-validateoutput,2)); 
        Denominator = sum(power(validatevalue(j,:)-ave_obs1,2));
        nse2 = 1 - Numerator/Denominator;

        [r21, rmse1] = rsquare(validateoutput,validatevalue(j,:));

        %figure
        %plot(validatevalue,validateoutput,'b.');
        %title(strcat(['R2 = ' num2str(r21) '; RMSE = ' num2str(rmse1)]))

        errors1=validatevalue(j,:)-validateoutput;
        trainSize1=size(validateoutput,2);
        for q=1:trainSize1
            zongshu2=sum(abs(validatevalue(j,q)));
            piancha2=sum(abs(validatevalue(j,q)-validateoutput(q)));
            pbias2(j,q)=100*piancha2/zongshu2;
        end
        pbi2=max(pbias2(j,:));
        %testoutput=mapminmax('reverse',y2,ts);
        %testvalue=mapminmax('reverse',testsample.t,ts);
        testoutput=y2;
        testvalue=testsamplet;

        ave_obs2 = sum(testvalue(:))/numel(testvalue);  
        Numerator1 = sum(power(testvalue-testoutput,2)); 
        Denominator1 = sum(power(testvalue-ave_obs2,2));
        nse3 = 1 - Numerator1/Denominator1;

        [r22, rmse2] = rsquare(testoutput,testvalue);
        trainSize2=size(testvalue,2);
        for q=1:trainSize2
            zongshu3=sum(abs(testvalue(q)));
            piancha3=sum(abs(testvalue(q)-testoutput(q)));
            pbias3(j,q)=100*piancha3/zongshu3;
        end
        pbi3=max(pbias3(j,:));

        %figure
        %plot(validatevalue,validateoutput,'b.');
        %title(strcat(['R2 = ' num2str(r21) '; RMSE = ' num2str(rmse1)]))

        erfang(j,1)=r2;
        erfang(j,2)=r21;
        erfang(j,3)=r22;
        erfang(j,4)=rmse;
        erfang(j,5)=rmse1;
        erfang(j,6)=rmse2;
        erfang(j,7)=nse1;
        erfang(j,8)=nse2;
        erfang(j,9)=nse3;
        erfang(j,10)=pbi1;
        erfang(j,11)=pbi2;
        erfang(j,12)=pbi3;


            if (r2>0.9&&r21>0.7)
                if rmse1<0.05
                break
                end
            end
        end
                
    zuheshu=zuheshu+1;
    %recording weight when the training and validation is ok
                    b1=net.b{1};
                    iw=net.iw{1,1};
                    b2=net.b{2};
                    lw=net.lw{2,1};
                    lw2=net.lw{3,2};
                    b3=net.b{3};
                    lw3=net.lw{4,3};
                    b4=net.b{4};
%                     lw4=net.lw{5,4};
%                     b5=net.b{5};
                    filename = strcat('weight-record',num2str(j),'.xlsx');
                    sheet = 1;
                    xlRange = 'A1:AX50';
                    xlswrite(filename,iw,sheet,xlRange);
                    xlRange = 'A52:A101';
                    xlswrite(filename,b1,sheet,xlRange);
                    xlRange = 'A104:AX138';
                    xlswrite(filename,lw,sheet,xlRange);
                    xlRange = 'A140:A174';
                    xlswrite(filename,b2,sheet,xlRange);
                    xlRange = 'A176:AX195';
                    xlswrite(filename,lw2,sheet,xlRange);
                    xlRange = 'A209:A228';
                    xlswrite(filename,b3,sheet,xlRange);
                    xlRange = 'A231:AX231';
                    xlswrite(filename,lw3,sheet,xlRange);
                    xlRange = 'A252:A252';
                    xlswrite(filename,b4,sheet,xlRange);
%                     xlRange = 'A281:T281';
%                     xlswrite(filename,lw4,sheet,xlRange);
%                     xlRange = 'A282:A282';
%                     xlswrite(filename,b5,sheet,xlRange);
                    xlRange = 'A283:O283';
                    xlswrite(filename,erfang(j,:),sheet,xlRange);

end
toc



