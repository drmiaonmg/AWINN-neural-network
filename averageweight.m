clc
tic
close all 
clear all 
rng default   

filename='weight-record-selected.xlsx';
sheet=2;
xlRange = 'Y1:AC1579';
Column11=xlsread(filename,sheet,xlRange); %查找输入
mun=[1579,868,756,610,433];
% filename='工作簿5.xlsx';
% sheet=1;
% xlRange = 'I1:I22';
% Column1=xlsread(filename,sheet,xlRange); %查找输入
% num=size(Column1,1);
% b11=zeros(50,num);
% w11=zeros(2100,num);
% b21=zeros(25,num);
% w21=zeros(1250,num);
% b31=zeros(1,num);
% w31=zeros(25,num);
for i=1:1
    huizong=zeros(mun(i),30*50);
        huizong2=zeros(mun(i),51*25);
        huizong3=zeros(mun(i),26*9);
        huizong4=zeros(mun(i),10*1);
    for j=1153:mun(i)
        
        a=Column11(j,i);
        filename= strcat('weight-record',num2str(a),'.xlsx');
        sheet=1;
        xlRange = 'A52:A101';
        b1=xlsread(filename,sheet,xlRange); %查找输入
        xlRange = 'A1:AC50';
        w1=xlsread(filename,sheet,xlRange); %查找输入
        w1size=size(w1);
        xlRange = 'A140:A164';
        b2=xlsread(filename,sheet,xlRange); %查找输入
        xlRange = 'A104:AX128';
        w2=xlsread(filename,sheet,xlRange); %查找输入
        w2size=size(w2);
        xlRange = 'A209:A218';
        b3=xlsread(filename,sheet,xlRange); %查找输入
        xlRange = 'A176:Y185';
        w3=xlsread(filename,sheet,xlRange); %查找输入
        w3size=size(w3);
        xlRange = 'A232:A232';
        b4=xlsread(filename,sheet,xlRange); %查找输入
        xlRange = 'A231:J231';
        w4=xlsread(filename,sheet,xlRange); %查找输入
        w4size=size(w4);
        for k=1:w1size(1)
            huizong(j,((30*(k-1)+1):(30*(k-1)+29)))=w1(k,:);
            huizong(j,30*k)=b1(k,1);
        end
        for k=1:w2size(1)     
            huizong2(j,((51*(k-1)+1):(51*(k-1)+50)))=w2(k,:);
            huizong2(j,51*k)=b2(k,1);
        end
        for k=1:w3size(1)     
            huizong3(j,((26*(k-1)+1):(26*(k-1)+25)))=w3(k,:);
            huizong3(j,26*k)=b3(k,1);
        end
        for k=1:w4size(1)     
            huizong4(j,((11*(k-1)+1):(11*(k-1)+10)))=w4;
            huizong4(j,11*k)=b4;
        end
    end
    
        
%     xlRange = 'A1:A50';
%     b1=xlsread(filename,sheet,xlRange); %查找输入
%     xlRange = 'A52:AH101';%50*34
%     w1=xlsread(filename,sheet,xlRange); %查找输入
%     w1size=size(w1);
%     xlRange = 'A104:A128';
%     b2=xlsread(filename,sheet,xlRange); %查找输入
%     xlRange = 'A156:AX180';%25*50
%     w2=xlsread(filename,sheet,xlRange); %查找输入
%     w2size=size(w2);
%     xlRange = 'A208:A208';
%     b3=xlsread(filename,sheet,xlRange); %查找输入
%     xlRange = 'A209:Y209';
%     w3=xlsread(filename,sheet,xlRange); %查找输入
%     w3size=size(w3);
%     b11(:,i)=b1;
%     for j=1:w1size(2)
%         w11((1+50*(j-1)):(50*j),i)=w1(:,j);
%     end
%     b21(:,i)=b2;
%     for j=1:w2size(2)
%         w21((1+25*(j-1)):(25*j),i)=w2(:,j);
%     end
%     b31(i)=b3;
%     for j=1:w3size(2)
%         w31(j,i)=w3(1,j);

    filename = 'Weight pooling.xlsx';
    sheet = 4*i;
    xlRange = 'A1:FZZ2000';
    xlswrite(filename,huizong,sheet,xlRange);
    sheet = 4*i-1;
    xlRange = 'A1:FZZ2000';
    xlswrite(filename,huizong2,sheet,xlRange);
    sheet = 4*i-2;
    xlRange = 'A1:FZZ2000';
    xlswrite(filename,huizong3,sheet,xlRange);
    sheet = 4*i-3;
    xlRange = 'A1:FZZ2000';
    xlswrite(filename,huizong4,sheet,xlRange);
    
    
end

    
    

