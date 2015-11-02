% FILENAME: tensor_facial_exp_eval1.m
% evaluate for arbitrary person expression by specifying 
% evaluate content x from unknown style and expression

clear
pack

load('TENSORFE_P10E3.mat')% 'c','u','s','rdim'
imagesize = 64*64;
CONTENT_PRECISION = 50;
centerNum = 8;

NEARZERO = 0.000000000001;
MINSQUARE = 1;
% % load test data
% load('CMU_FE_IMAGEp4e3'); % FEIMAGE{pid,expid}
%                             % expression_pid
%                             % expression_index

% arbitrary test data
expression_pid_t = 6;
expression_index_t = [26 34;54 64;67 75];


load('CMUFEIMAGE13'); % CMUFEIMAGE

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
id = 1; % C
expression_index{id} = [53 64; 20 30; 10  18]; % 12 11 9
expression_pid{id} = 3;
id = 2; % D
expression_index{id} = [19 33; 52  66; 68 75]; % 15 15 8
expression_pid{id} = 4;
id = 3; % E
expression_index{id} = [19 33 ; 55 64; 67 75]; % 15 10 9
expression_pid{id} = 5;
id = 4; % K
expression_index{id} = [64 75; 7 19; 23 47]; % 12 13 25
expression_pid{id} = 11;
id = 5; % A
expression_index{id} = [19 28; 33 43; 46 63]; % 10 11 18
expression_pid{id} = 1;
id = 6; % F
expression_index{id} = [25 33; 55 64; 67 75]; % 9 10 9
expression_pid{id} = 6;
id = 7; % H
expression_index{id} = [22 33; 66 75; 46 59]; % 12 10 14
expression_pid{id} = 8;
id = id+1; % J
expression_index{id} = [60 73; 34 53; 14 29]; % 14 20 16
expression_pid{id} = 10;
id = id+1; % L
expression_index{id} = [54 64; 7  16; 23 35]; % 11 10 13
expression_pid{id} = 12;
id = id+1; % M
expression_index{id} = [55 63; 7  16; 26 34]; % 9 10 9
expression_pid{id} = 13;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FEImage = double(CMUFEIMAGE)/256;
ALLFRAME_NUM = 75;
EXP_NUM = 3;

numofstyle = size(c,1);
numofexp = size(c,2);
pixelDIM = size(c,3);

id_num = id;
exp_num = EXP_NUM;


Tp =tmul(c, u{3}',3);

a_ds = 5;
b_ds = 5;

iterationNum= 7;

style_mean = u{1};
exp_mean = u{2};
        
THRESHOLD_S = 0.8;
THRESHOLD_E = 0.7;

% for i = 1: numofstyle

% for i = 1: id_num
for i = 1: 1
%     for j = 1: numofexp
%     for j = 1:3       
        sindex = 1:numofstyle;
        ws(sindex,1) = 1/length(sindex);
        eindex = 1:numofexp;
        we(eindex,1) = 1/length(eindex);
        
        pid = expression_pid{i};
        TestData = FEImage((pid-1)*ALLFRAME_NUM +1:(pid)*ALLFRAME_NUM,:)';
%         TestData = FEIMAGE{i}{j}';
        
        cycleframenum = size(TestData,2);
        ci = [1:cycleframenum];
        
%         % true content
%         x = ci*2*pi/cycleframenum;
%         xt = [cos(x') sin(x')];
        

        xt = zeros(ci,2);
        
        
        %%%%%%%%%%%%%%%%%%%%%
        for j=1:3
            testcycleframenum = (expression_index{i}(j,2) - expression_index{i}(j,1)) +1;
            ti = 1:testcycleframenum;
			xx=ti*2*pi/testcycleframenum;
            
            eval1=[cos(xx') sin(xx')];         
            
            xt(expression_index{i}(j,1):expression_index{i}(j,2),:) = eval1;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%
        
        
        
        xef = [];
        xef_error = [];
        

        styleindexe = [];
        stylevaluee = [];
        expindexe = [];
        expvaluee = [];
        
        for k = 1:cycleframenum

            xe = [];
            xe_error = [];
            xe_new = [];
            xe_new_error = [];        

            a_sigma = 50000;
            b_sigma = 50000;
        
            wse = [];
            wee = [];
            
            avgse = [];
            avgee = [];
            
            avgse_e = [];
            avgee_e = []; 
        
%             imagev = FEIMAGE{i}{j}(k,:)';
            imagev = TestData(:,k);
            
            sindex = 1:numofstyle;
            ws(sindex,1) = 1/length(sindex);
            eindex = 1:numofexp;
            we(eindex,1) = 1/length(eindex);
        
        
            for iter =1:iterationNum
                sc = ws'*style_mean;
                
                ec = we'*exp_mean;
	
              
                Tpsc = tmul(Tp, sc',1);
                Tpscec = tmul(Tpsc, ec',2);           
                CFe = squeeze(Tpscec);
                
                rCFe = reshape(CFe, (centerNum+3), imagesize);
            
            
                [evalx, errorx] = findcontent_fullsearch(imagev, rCFe, centerNum, CONTENT_PRECISION);
                
                xe = [xe; evalx];
                
%                 xe_er = (xt(k,:)-evalx);
%                 xe_er2 = sum(xe_er.*xe_er);
%                 
%                 xe_error = [xe_error; xe_er2];  
%                           
                
                %%% compute error in each view model %%%%
                expFactor = [];
                expError = [];
                for jj = 1:numofexp
                    
     			    Tpscvi = tmul(Tpsc,exp_mean(jj,:)',2);
                    
                    newCFi = squeeze(Tpscvi);
                    
                    % computation of similarity
                    rnewCFi = reshape(newCFi, (centerNum+3), imagesize);
                    
                    [error_distance, factor] = measure_reconst_one_image(imagev, rnewCFi,centerNum, evalx, a_sigma);
                    
					expFactor = [expFactor; factor];
                    expError = [expError; error_distance];
                end % for j = 1
                        
                % simple regression model based on estimated probability ratio
                sumweight = sum(expFactor);
                
                if sumweight >=NEARZERO
                    we = expFactor/sumweight;
                end
                avge = we'*expError;
	
                ec = we'*exp_mean;            
                
                %% error estimation way2
                CF_newexp = squeeze(tmul(Tpsc,ec',2));
                rCF_newexp = reshape(CF_newexp, (centerNum+3), imagesize);
                    
                [avge_e] = compareimage_content(imagev, rCF_newexp,centerNum, evalx);            
                
                        
                Tpve = tmul(Tp, ec',2);
                Tpvesc = tmul(Tpve, sc',1);
                
              
                CF_new = squeeze(Tpvesc);
                nCF_new = reshape(CF_new, (centerNum+3), imagesize);
                
                evalx_new = findcontent_fullsearch(imagev, nCF_new,centerNum, CONTENT_PRECISION);
                
                xe_new = [xe_new; evalx_new];
                
%                 xe_new_er = (xt(k,:)-evalx_new);
%                 xe_new_er2 = sum(xe_new_er.*xe_new_er);
%                 
%                 xe_new_error = [xe_new_error; xe_new_er2];  
                
                %%%%%%%% re-estimation for style too %%%%%%%% 
                styleFactor = [];
                styleError = [];
                for ii = 1:numofstyle
                            
                    newCFe = squeeze(tmul(Tpve,style_mean(ii,:)',1));
                    
                    % computation of similarity
                    rnewCFe = reshape(newCFe, (centerNum+3), imagesize);
                    
                    [error_distance2, factor2] = measure_reconst_one_image(imagev, rnewCFe,centerNum, evalx_new, b_sigma);
                    
					styleFactor = [styleFactor; factor2];
                    styleError = [styleError; error_distance2];
                    
                end % for ii = 1
                
                % simple regression model based on estimated probability ratio
                sumweight = sum(styleFactor);
                if sumweight >=NEARZERO
                    ws = styleFactor/sumweight;
                end
                
                avgs = ws'*styleError;
	
                sc = ws'*style_mean;
                
                CF_newstyle = squeeze(tmul(Tpve,sc',1));
                rCF_newstyle = reshape(CF_newstyle, (centerNum+3), imagesize);
                    
                [avgs_e] = compareimage_content(imagev, rCF_newstyle,centerNum, evalx_new);	
	
                        
                wee = [wee we];
                wse = [wse ws];
                
                avgse = [avgse avgs];
                avgee = [avgee avge];                    
                
                avgse_e = [avgse_e avgs_e];
                avgee_e = [avgee_e avge_e];                       
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if a_sigma>MINSQUARE
                    a_sigma = a_sigma /a_ds;
                end
               
                if b_sigma>MINSQUARE
                    b_sigma = b_sigma /b_ds;
                end

            end %iter

            [estylevalue, estyleindex] = max(ws);
            
            if estylevalue > THRESHOLD_S
                styleindexe = [styleindexe; estyleindex];
            else
                styleindexe = [styleindexe; 0];
%                 evalx_new = [0 0];
            end
            
            stylevaluee = [stylevaluee; estylevalue];
            
            [eexpvalue, eexpindex] = max(we);
            
            if eexpvalue > THRESHOLD_E
                expindexe = [expindexe; eexpindex];
            else
                expindexe = [expindexe; 0];
                evalx_new = [0 0];
            end
            
            expvaluee = [expvaluee; eexpvalue];            
            
            xef = [xef; evalx_new];
%             xef_error = [xef_error;sqrt(xe_new_er2)];
            
%             figure
%             subplot(1,4,1)
%             plot(avgee,'r.-')
%             hold on
%             plot(avgse,'b.-')
% 	
%             subplot(1,4,2)
%             plot(avgee_e,'r*-')
%             hold on
%             plot(avgse_e,'b*-')
%             
%             subplot(1,4,3)
%             plot(wse')
%             
%             subplot(1,4,4)
%             plot(wee')
% 	
%             
%             
%             tss = strcat('P: ',int2str(i),' F=',int2str(k));
% %             tss = strcat('P: ',int2str(i),' E:',int2str(j),' F=',int2str(k));
%             title(tss)
%                            
%             drawnow 

        end
        
            
        figure
        subplot(1,2,1)
		plot(xt(:,1),'bo-')
        hold on
		plot(xef(:,1),'ro-')
        legend('true x value','estimated x value')        
        xlabel('frame number')
        
        subplot(1,2,2)
		plot(xt(:,2),'b*-')
        hold on
		plot(xef(:,2),'r*-')

% 		title('correct and estimated content vector')
		legend('true y value','estimated y value')
		xlabel('frame number')
	% 			plot(xtd(:,1),'go-')
	% 			plot(xtd(:,2),'g*-')
	% 			legend('est. x','est.y','al. x','al. y','Nal x','Nal. y')
    
    
    
        figure
        intensity_true = (1-xt(:,2))/2;
        intensity_estimated = (1-xef(:,2))/2;
        
        plot(intensity_true, 'b.-')
        hold on
        plot(intensity_estimated, 'r*-')
        
        figure       
        plot(styleindexe,'r.-')
        hold on
        plot(expindexe,'b*-')        
        
        
%         bar(xef_error)   
% %         ts = strcat('person id =', int2str(i), ' expression id = ',int2str(j))
        ts = strcat('person id =', int2str(i))
        title(ts)
%     end
		figure
		plot(expvaluee, 'r.-')
		grid on
		hold on
		plot(stylevaluee,'bo-')
		legend('exp prob','style prob')
        title('Prob. for maximum expression and style class')
end
    
