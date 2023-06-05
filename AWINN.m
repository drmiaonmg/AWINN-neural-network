clc
 
clear all
 
close all

filename='data.xlsx';
sheet=1;
xlRange = 'A2:AC67';
Column1=xlsread(filename,sheet,xlRange); %查找输入
filename='data.xlsx';
sheet=1;
xlRange = 'A2:M10001';
Column1=xlsread(filename,sheet,xlRange); %查找输入
filename='data.xlsx';
sheet=1;
xlRange = 'AD2:AF67';
Column2=xlsread(filename,sheet,xlRange); %查找输入

% 
%  
% 
p=Column1';
%  
t=Column2(:,1)';

 
[p1,ps]=mapminmax(p,0,1);
 
[t1,ts]=mapminmax(t,0,1);
%
%[trainInd,valInd,testInd] = dividerand(Q,trainRatio,valRatio,testRatio)
% i=1:5;
% num=3;
% y=combntns(i,num);
% trainsample.p=p1(:,1:26);
% valsample.p=p1(:,27:31);
testsamplep=p1(:,61:66);
% trainsample.t=t1(:,1:26);
%valsample.t=t1(:,27:31);
testsamplet=t1(:,61:66);

func={'trainlm','traingdx','traingda','traingd','trainbr','trainscg'};
func1=size(func);
trans={'tansig','logsig'};
trans1=size(trans);
% for hiddenLayerSize = 12:18
%     for o=1:6
%         trainFcn = char(func(o));  % Bayesian Regularization
%         for p=1:1
%             net.layers{1}.transferFcn = 'logsig';
%             for j=1:10
%                 c=setdiff(i,y(j,:));
%                for k=1:3
%                     trainsample.p=[trainsample.p p1(:,((y(j,k)-1)*101+1):(y(j,k)*101))];
%                     trainsample.t=[trainsample.t t1(:,((y(j,k)-1)*101+1):(y(j,k)*101))];
%  %               end
%                 for l=1:2 
%                     if l==1
%                         valsample.p=p1(:,((c(l)-1)*101+1):c(l)*101);
%                         testsample.p=p1(:,((c(l+1)-1)*101+1):c(l+1)*101);
%                         valsample.t=t1(:,(c(l)-1)*101+1:c(l)*101);
%                         testsample.t=t1(:,((c(l+1)-1)*101+1):c(l+1)*101);
%                     else
%                         valsample.p=p1(:,((c(l)-1)*101+1):c(l)*101);
%                         testsample.p=p1(:,((c(l-1)-1)*101+1):c(l-1)*101);
%                         valsample.t=t1(:,((c(l)-1)*101+1):c(l)*101);
%                         testsample.t=t1(:,((c(l-1)-1)*101+1):c(l-1)*101);
%                     end
                %elseif i==2

                        %trainsample.p=p1(:,206:505);
                    %valsample.p=p1(:,1:100);
                    %testsample.p=p1(:,401:505);
                    %trainsample.t=t1(:,101:400);
                    %valsample.t=t1(:,1:100);
                    %testsample.t=t1(:,401:505);
x = p1;
p2=p1(:,31:60);
t = t1;
t2=t1(:,31:60);
fenzu=[11 2 3 9 4 ];
%zuheshu=fenzu(1)*fenzu(2)*fenzu(3)*fenzu(4)*fenzu(5);
zuheshu=0;
% zup1=1:11;
% 
% zup2=12:13;
% 
% zup3=14:16;
% 
% zup4=17:25;
% 
% zup5=26:29;
n=1:30;
% for ju1=1:11
%     zuhe=[ju1 30];
%     for ju2=12:13
%         zuhe=[zuhe ju2];
%         for ju3=14:16
%             zuhe=[zuhe ju3];
%             for ju4=17:25
%                 zuhe=[zuhe ju4];
%                 for ju5=26:29
%                     zuhe=[zuhe ju5];
% 
%                     zuhe=[zuhe 30];
                    % m=1:5;

%                     num=6;
%                     y=combntns(6,num)
                    m=1:5;
                    num=4;
                    y=combntns(m,num);
                    % k=3;
                    % trainxulie=[((y(k,1)-1)*6+1):y(k,1)*6 ((y(k,2)-1)*6+1):y(k,2)*6 ((y(k,3)-1)*6+1):y(k,3)*6 ((y(k,4)-1)*6+1):y(k,4)*6]
%                         trainxulie=zeros(5,24);
%                         valxulie=zeros(5,6);
%                         trainvaluep=zeros(145,24);
%                         trainvalue=zeros(5,24);
%                         valvaluep=zeros(145,6);
%                         validatevalue=zeros(5,6);
%                         trainSize=zeros(1,5);
%                         trainSize1=zeros(1,5);
%                         trainSize2=zeros(1,5);
%                         pbias1=zeros(5,24);
%                         pbias2=zeros(5,6);
%                         pbias3=zeros(5,6);
%                         erfang=zeros(5,12);
% % xlgoal=[0.002,0.015,0.027,0.015,0.01,0.001,0.01,0.01,0.04,0.01;0.002,0.025,0.035,0.0335,0.0335,0.001,0.032,0.0335,0.03,0.0335;0.003,0.0325,0.02,0.031,0.03,0.001,0.02,0.005,0.0112,0.025;0.0015,0.035,0.0325,0.04,0.04,0.001,0.02,0.0375,0.0375,0.0375];
% xlgoal=[0.001,0.015,0.015,0.015,0.015,0.001,0.018,0.008,0.015,0.02;0.001,0.03,0.02,0.01,0.015,0.001,0.015,0.02,0.028,0.015;0.0015,0.028,0.035,0.01,0.015,0.001,0.018,0.03,0.0112,0.015;0.0015,0.015,0.015,0.015,0.015,0.001,0.015,0.015,0.015,0.02];
% for i=1:4
% filename = 'nh4jieguo.xlsx';
% sheet = i;
% xlRange ='A1:AC20';
% huizong=xlsread(filename,sheet,xlRange);
% % sheet = 3;
% xlRange ='A104:T116';
% huizong1=xlsread(filename,sheet,xlRange);
% % sheet = 5;
% xlRange ='A176:M181';
% huizong2=xlsread(filename,sheet,xlRange);
% % sheet = 7;
% xlRange = 'A231:F231';
% huizong3=xlsread(filename,sheet,xlRange);
% % sheet = 2;
% xlRange = 'A52:A71';
% huizongb=xlsread(filename,sheet,xlRange);
% % sheet = 4;
% xlRange =  'A140:A152';
% huizongb1=xlsread(filename,sheet,xlRange);
% % sheet = 6;
% xlRange = 'A209:A214';
% huizongb2=xlsread(filename,sheet,xlRange);
% % sheet = 8;
% xlRange = 'A252:A252';
% huizongb3=xlsread(filename,sheet,xlRange);
for k=1:5
    trainxulie(k,1:24)=[((y(k,1)-1)*6+1):y(k,1)*6 ((y(k,2)-1)*6+1):y(k,2)*6 ((y(k,3)-1)*6+1):y(k,3)*6 ((y(k,4)-1)*6+1):y(k,4)*6 ];
    valxulie(k,1:6)=setdiff(n,trainxulie(k,1:24));
    trainvaluep(((k-1)*29+1):k*29,1:24)=p1(1:29,trainxulie(k,1:24));
    valvaluep(((k-1)*29+1):k*29,1:6)=p1(1:29,valxulie(k,1:6));
    trainvalue(k,1:24)=t1(1,trainxulie(k,1:24));
    validatevalue(k,1:6)=t1(1,valxulie(k,1:6));
end
for k=6:10
    trainxulie(k,1:24)=[((y(k-5,1)-1)*6+1):y(k-5,1)*6 ((y(k-5,2)-1)*6+1):y(k-5,2)*6 ((y(k-5,3)-1)*6+1):y(k-5,3)*6 ((y(k-5,4)-1)*6+1):y(k-5,4)*6 ];
    valxulie(k,1:6)=setdiff(n,trainxulie(k,1:24));
    trainvaluep(((k-1)*29+1):k*29,1:24)=p2(1:29,trainxulie(k,1:24));
    valvaluep(((k-1)*29+1):k*29,1:6)=p2(1:29,valxulie(k,1:6));
    trainvalue(k,1:24)=t2(1,trainxulie(k,1:24));
    validatevalue(k,1:6)=t2(1,valxulie(k,1:6));
end
%select different dataset to retrain models and set different values of MSE
xlgoal=[0.0063,0.001,0.0055,0.0055,0.007,0.0063,0.006,0.0055,0.0055,0.007];
for i=1:1
filename = 'Weight pooling.xlsx';
sheet = 8*i;
xlRange ='A1:AC50';
huizong=xlsread(filename,sheet,xlRange);
sheet = 8*i-1;
xlRange ='A104:AX128';
huizong1=xlsread(filename,sheet,xlRange);
sheet = 8*i-2;
xlRange ='A176:Y185';
huizong2=xlsread(filename,sheet,xlRange);
sheet = 8*i-3;
xlRange = 'A231:J231';
huizong3=xlsread(filename,sheet,xlRange);
sheet = 8*i-4;
xlRange = 'A52:A101';
huizongb=xlsread(filename,sheet,xlRange);
sheet = 8*i-5;
xlRange =  'A140:A164';
huizongb1=xlsread(filename,sheet,xlRange);
sheet = 8*i-6;
xlRange = 'A209:A218';
huizongb2=xlsread(filename,sheet,xlRange);
sheet = 8*i-7;
xlRange = 'A232:A232';
huizongb3=xlsread(filename,sheet,xlRange);

% trainvaluep(1:29,1:30)=p1(1:29,1:30);
% trainvalue(1,1:30)=t1(1,1:30);
erfang=zeros(6,12);
                    for k=1:10  %选取不同的数据集
%                         trainxulie=zeros(5,24);
%                         valxulie=zeros(5,6);
%                         trainvaluep=zeros(145,24);
%                         trainvalue=zeros(5,24);
%                         valvaluep=zeros(145,6);
%                         validatevalue=zeros(5,6);
%                         trainSize=zeros(1,5);
%                         net.b{1}=huizongb;
%                         net.iw{1,1}=huizong;
%                         net.b{2}=huizongb1;
%                         net.lw{2,1}=huizong1;
%                         net.lw{3,2}=huizong2;
%                         net.b{3}=huizongb2;
%                        net.lw{4,3}=huizong3;
%                         net.b{4}=huizongb3;
%                         net.lw{5,4};
%                         b5=net.b{5};
%                         trainSize1=zeros(1,5);
%                         trainSize2=zeros(1,5);
                        pbias1=zeros(6,24);
                        pbias2=zeros(6,6);
                        pbias3=zeros(6,6);
                        
                        pbi1=zeros(6,1);
                        pbi3=zeros(6,1);
                        r2=zeros(6,1);
                        rmse=zeros(6,1);
                        nse=zeros(6,1);
                        r21=zeros(6,1);
                        rmse1=zeros(6,1);
                        nse1=zeros(6,1);                        
                        r22=zeros(6,1);
                        rmse2=zeros(6,1);
                        nse3=zeros(6,1);
                        
                        
%                  if k==1%trainxulie(k,1:24)
%                     trainxulie1=[((y(k,1)-1)*6+1):y(k,1)*6 ((y(k,2)-1)*6+1):y(k,2)*6 ((y(k,3)-1)*6+1):y(k,3)*6 ((y(k,4)-1)*6+1):y(k,4)*6];
%                     valxulie1=setdiff(n,trainxulie1);
%                     trainvaluep(((k-1)*29+1):k*29,1:24)=p1(1:29,trainxulie1);
%                     valvaluep(((k-1)*29+1):k*29,1:6)=p2(1:29,valxulie1);
%                     trainvalue(k,1:24)=t1(1,trainxulie1);
%                     validatevalue(k,1:6)=t2(1,valxulie1);
% 
%                     trainxulie1=trainxulie(k,:);
%                     valxulie1=valxulie(k,:);
                   for j=1:1
                    trainFcn = 'trainscg';  % Bayesian Regularization
                    hiddenLayerSize =50;
                    hiddenLayerSize1=25;
                    hiddenLayerSize2=10;
                    hiddenLayerSize3=10;
                    %net = feedforwardnet (hiddenLayerSize,trainFcn);
                    net = feedforwardnet ([hiddenLayerSize,hiddenLayerSize1,hiddenLayerSize2],trainFcn);
                    % Setup Division of Data for Training, Validation, Testing
                    %RandStream.setGlobalStream(RandStream('mt19937ar','seed',1)); % to get constant result
                    net.divideFcn = 'divideind'; % Divide targets into three sets using blocks of indices
                    %[trainInd,valInd,testInd] = divideind(36,1:30,0:0,31:36);%需要标记k
                    trainInd=1:24;
%                     valInd= valxulie(k,1:6);
%                     testInd=31:36;
                    % net.divideParam.trainInd = 1;
                    % net.divideParam.valRatio = 0;
                    % net.divideParam.testRatio = 0;
                    net.trainParam.showwindow=0;
                    net.divideParam.trainInd=trainInd;
%                     net.divideParam.valInd=valInd;
%                     net.divideParam.testInd=testInd;
                    %net.numInputs=29;
                    net=configure(net,trainvaluep(((k-1)*29+1):k*29,1:24),trainvalue(k,1:24));
                    net.b{1}=huizongb;
                    appldd=net.b{1};
                    net.iw{1,1}=huizong;
                    net.b{2}=huizongb1;
                    net.lw{2,1}=huizong1;
                    net.lw{3,2}=huizong2;
                    net.b{3}=huizongb2;
                    net.lw{4,3}=huizong3;
                    net.b{4}=huizongb3;
                    %TRAINING PARAMETERS
                    net.trainParam.show=1;  % of ephocs in display
                    net.trainParam.lr=0.005;  % learning rate
                    net.trainParam.max_fail=800;
                    net.trainParam.epochs=500000;  % max epochs
                    net.trainParam.goal=xlgoal(k);  % training goal
%                     net.trainParam.goal=0.0075;  % training goal
                    net.performFcn='mse';  % Name of a network performance function %type help nnperformance
                    %net.numLayers = 3;  % 创建更多的隐藏层以及添加激活层函数
                    %net.layerConnect(3,1) = 1;
                    %net.outputConnect = [0 0 1 ];
                    %net.layers{1}.size = 12; %第二层有八个节点
                    net.layers{1}.transferFcn = 'logsig';
                    %net.layers{2}.size = 20; %第二层有八个节点
                    net.layers{2}.transferFcn = 'logsig';
                     net.layers{3}.transferFcn = 'logsig';
                     %net.layers{3}.transferFcn = 'logsig';
                    %net.outputs{2}.range = [1 20];
                    % Train the Network
                    %循环

                    [net,tr] = train(net,trainvaluep(((k-1)*29+1):k*29,1:24),trainvalue(k,1:24)); 
                    % Test the Network
                    y0 = net(trainvaluep(((k-1)*29+1):k*29,1:24));
                    e = gsubtract(trainvalue(k,1:24),y0);
                    performance = perform(net,trainvalue,y0);
                    y1 = net(valvaluep(((k-1)*29+1):k*29,:));
                    performance1 = perform(net,validatevalue,y1);
                    y2 = net(testsamplep);
                    performance2 = perform(net,testsamplet,y2);
                    % View the Network

                    trainoutput=mapminmax('reverse',y0,ts);
                    trainvalue1=mapminmax('reverse',trainvalue(k,:),ts);%正常的验证数据
                    %trainoutput=y0;
                    %trainvalue=trainsample.t;%正常的验证数据
        
                    ave_obs = sum(trainvalue1)/numel(trainvalue1);   %实测数据平均数
                    Numerator = sum(power(trainvalue1-trainoutput,2));  %分子
                    Denominator = sum(power(trainvalue1-ave_obs,2));%分母
                    nse(k) = 1 - Numerator/Denominator;%实测数据纳什系数

                    [r2(k), rmse(k)] = rsquare(trainoutput,trainvalue1);
                    errors=trainvalue1-trainoutput;
                    trainSize=size(trainoutput,2);
                    for q=1:trainSize
                        zongshu1=abs(trainvalue1(q));
                        piancha1=abs(trainvalue1(q)-trainoutput(q));
                        pbias1(k,q)=100*piancha1/zongshu1;
                    end
                    pbi1(k)=max(pbias1(k,:));
%                     filename = strcat('pianchajilu',num2str(i),'.xlsx');
%                     sheet = 1;
%                     xlRange =strcat('A',num2str(k),':X',num2str(k));
%                     xlswrite(filename,pbias1(k,:),sheet,xlRange);
                   end
%                     validateoutput=y1;
%                     validatevalue=valsample.t;%正常的验证数据
                    validateoutput=mapminmax('reverse',y1,ts);
                    validatevalue1=mapminmax('reverse',validatevalue(k,:),ts);%正常的验证数据

                    ave_obs1 = sum(validatevalue1)/numel(validatevalue1);   %实测数据平均数
                    Numerator = sum(power(validatevalue1-validateoutput,2));  %分子
                    Denominator = sum(power(validatevalue1-ave_obs1,2));%分母
                    nse2(k) = 1 - Numerator/Denominator;%实测数据纳什系数

                    [r21(k), rmse1(k)] = rsquare(validateoutput,validatevalue1);

                    %figure
                    %plot(validatevalue,validateoutput,'b.');
                    %title(strcat(['R2 = ' num2str(r21) '; RMSE = ' num2str(rmse1)]))

                    errors1=validatevalue1-validateoutput;
                    trainSize1=size(validateoutput,2);
                    for q=1:trainSize1
                        zongshu2=abs(validatevalue1(q));
                        piancha2=abs(validatevalue1(q)-validateoutput(q));
                        pbias2(k,q)=100*piancha2/zongshu2;
                    end
                    pbi2(k)=max(pbias2(k,:));
%                     filename = strcat('pianchajilu',num2str(i),'.xlsx');
%                     sheet = 2;
%                     xlRange =strcat('A',num2str(k),':F',num2str(k));
%                     xlswrite(filename,pbias2(k,:),sheet,xlRange);
                    testoutput=mapminmax('reverse',y2,ts);
                    testvalue=mapminmax('reverse',testsamplet,ts);%正常的验证数据
%                     testoutput=y2;
%                     testvalue=testsamplet;%正常的验证数据

                    ave_obs2 = sum(testvalue)/numel(testvalue);   %实测数据平均数
                    Numerator1 = sum(power(testvalue-testoutput,2));  %分子
                    Denominator1 = sum(power(testvalue-ave_obs2,2));%分母
                    nse3(k) = 1 - Numerator1/Denominator1;%实测数据纳什系数

                    [r22(k), rmse2(k)] = rsquare(testoutput,testvalue);
                    trainSize2=size(testvalue,2);
                    for q=1:trainSize2
                        zongshu3=sum(abs(testvalue(q)));
                        piancha3=sum(abs(testvalue(q)-testoutput(q)));
                        pbias3(k,q)=100*piancha3/zongshu3;
                    end
                    pbi3(k)=max(pbias3(k,:));
                    filename = strcat('pianchajilucod',num2str(i),'.xlsx');
                    sheet = 3;
                    xlRange =strcat('A',num2str(k),':F',num2str(k));
                    xlswrite(filename,testoutput,sheet,xlRange);
                    %figure
                    %plot(validatevalue,validateoutput,'b.');
                    %title(strcat(['R2 = ' num2str(r21) '; RMSE = ' num2str(rmse1)]))

                    erfang(k,1)=r2(k);
                    erfang(k,2)=r21(k);
                    erfang(k,3)=r22(k);
                    erfang(k,4)=rmse(k);
                    erfang(k,5)=rmse1(k);
                    erfang(k,6)=rmse2(k);
                    erfang(k,7)=nse(k);
                    erfang(k,8)=nse2(k);
                    erfang(k,9)=nse3(k);
                    erfang(k,10)=pbi1(k);
                    erfang(k,11)=pbi2(k);
                    erfang(k,12)=pbi3(k);
                    
                    
%                         if (r2>0.9&&r21>0.7)&&r22>0.7
%                             if rmse1<0.2&&rmse2<0.2
%                             break
%                             end
%                                 
%                         end
                    
                    zuheshu=zuheshu+1;
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
%record new weights
                    filename = strcat('cod',num2str(i),'.xlsx');
                    sheet = k;
                    xlRange ='A1:AC50';
                    xlswrite(filename,iw,sheet,xlRange);
                    xlRange = 'A52:A101';
                    xlswrite(filename,b1,sheet,xlRange);
                    xlRange ='A104:AX128';
                    xlswrite(filename,lw,sheet,xlRange);
                    xlRange ='A140:A164';
                    xlswrite(filename,b2,sheet,xlRange);
                    xlRange = 'A176:Y185';
                    xlswrite(filename,lw2,sheet,xlRange);
                    xlRange =  'A209:A218';
                    xlswrite(filename,b3,sheet,xlRange);
                    xlRange = 'A231:J231';
                    xlswrite(filename,lw3,sheet,xlRange);
                    xlRange = 'A232:A232';
                    xlswrite(filename,b4,sheet,xlRange);
%                     xlRange = 'A281:T281';
%                     xlswrite(filename,lw4,sheet,xlRange);
%                     xlRange = 'A282:A282';
%                     xlswrite(filename,b5,sheet,xlRange);
%                     xlRange = 'A283:O289';
%                     xlswrite(filename,erfang,sheet,xlRange);
%                     elseif k==2
%                         
%                     elseif k==3
%                         
%                     elseif k==4
%                         
%                     else
%                  end
                    
%                     filename = strcat('zuihou333chonhuansilu',num2str(i),'.xlsx');
%                     sheet = k;
%                     xlRange = 'A283:O289';
%                     xlswrite(filename,erfang,sheet,xlRange);
                    end
filename = strcat('coderfang',num2str(6),'scg.xlsx');
                    sheet = i;
                    xlRange = 'A283:O292';
                    xlswrite(filename,erfang,sheet,xlRange);
end
