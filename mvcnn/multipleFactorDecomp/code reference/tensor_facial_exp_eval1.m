% FILENAME: tensor_facial_exp_eval1.m
% evaluate content x from unknown style and expression

clear
pack

load('TENSORFE_P4E3.mat')% 'c','u','s','rdim'
imagesize = 64*64;
CONTENT_PRECISION = 50;
centerNum = 9;

% load test data
load('CMU_FE_IMAGEp4e3')% FEIMAGE{pid,expid}

numofstyle = size(c,1);
numofexp = size(c,2);
pixelDIM = size(c,3);

Tp =tmul(c, u{3}',3);

a_ds = 3;
b_ds = 3;

iterationNum= 7;

style_mean = u{1};
exp_mean = u{2};
        
for i = 1: numofstyle
    for j = 1: numofexp
% for i = 1: 1
%     for j = 1: numofexp        
        sindex = 1:numofstyle;
        ws(sindex,1) = 1/length(sindex);
        eindex = 1:numofexp;
        we(eindex,1) = 1/length(eindex);
        
        TestData = FEIMAGE{i}{j}';
        
        cycleframenum = size(TestData,2);
        ci = [1:cycleframenum];
        
        % true content
        x = ci*2*pi/cycleframenum;
        xt = [cos(x') sin(x')];
        
        xef = [];
        xef_error = [];
        

        
        xe = [];
        xe_error = [];
        xe_new = [];
        xe_new_error = [];        

        a_sigma = 1000;
        b_sigma = 1000;
        

       
        wse = [];
        wee = [];
        
        avgse = [];
        avgee = [];

        for k = 1:cycleframenum

            imagev = FEIMAGE{i}{j}(k,:);
            
            sc = ws'*style_mean;
            
            ec = we'*exp_mean;

%             Tpsi = tmul(Tp, style_mean(i,:)',1);
            Tpsc = tmul(Tp, sc',1);
%             Tpsiec = tmul(Tpsi, ec',2);
            Tpscei = tmul(Tpsc, exp_mean(j,:)',2);            
            CFe = squeeze(Tpscei);
            
            rCFe = reshape(CFe, (centerNum+3), imagesize);
        
        
            [evalx, errorx] = findcontent_fullsearch(imagev', rCFe, centerNum, CONTENT_PRECISION);
            
            xe = [xe; evalx];
            
            xe_er = (xt(k,:)-evalx);
            xe_er2 = sum(xe_er.*xe_er);
            
            xe_error = [xe_error; xe_er2];  
                      
            
        end
        figure
        subplot(1,2,1)
		plot(xt(:,1),'bo-')
        hold on
		plot(xt(:,2),'b*-')
		plot(xe(:,1),'ro-')
		plot(xe(:,2),'r*-')

		title('correct and estimated content vector')
		legend('true x value','true y value','estimated x value','estimated y value')
		xlabel('frame number')
	% 			plot(xtd(:,1),'go-')
	% 			plot(xtd(:,2),'g*-')
	% 			legend('est. x','est.y','al. x','al. y','Nal x','Nal. y')
        subplot(1,2,2)
        bar(xe_error)    
        ts = strcat('person id =', int2str(i), ' expression id = ',int2str(j))
        title(ts)
    end

end
    