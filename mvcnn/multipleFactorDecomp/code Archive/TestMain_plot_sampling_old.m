function TestMain_plot(p_data_file,p_use_sampling,p_use_wts)
% FILENAME: tensor_facial_exp_eval1.m
% This experimnet for speaker identification
%train over 26x2x10, test for 26x1x10.
if nargin<3
    p_use_wts=0;
end
if nargin<2
    p_use_sampling=0;
end
%clear
pack
%diary('run_log_HoG-3-52-50.txt');
%diary('run_log_Gabor-set2_26-20-(1-2).txt');
diary 'run_log_test'
echo on
%load('TEST_DATA.mat')% 'c','u','s','rdim'
disp 'train over reptns 1,2 and test over reptn 3'
subjects = {'*'}
letters = {'*'}
repeats = {'3'}
[tSeqs,tSizes,tLetters,tSubjects] = LoadDataAVLetterFiles(letters,subjects,repeats);


if nargin<1
    % % load test data
    load('AVLetter_Learn_HOSVD_52-10-(1-2)');
else
    load(p_data_file);
end
%c core tensor
%u{1} letter bases
%u{2} person style bases
threshold = 1E-10;
numofletters = size(tSeqs,1);%size(c,1);
numofstyles = size(tSeqs,3);%size(c,2);
pixelDIM = size(c,3);
letter_bases = u{1};
style_bases = u{2};

ltr_wts = 2*(length(letter_bases):1);
ltr_wts = ltr_wts/sum(ltr_wts);
%Tp =tmul(c, u{3}',3);
style_mean = mean(style_bases,1)';% estimate the initial of person style vectors from the trained style basis vectors 
letter_mean = mean(letter_bases,1)';
ltr_labels = [1:26,1:26]';
style_labels = (1:10)';
fig = plotting(letter_bases(:,1:3),labels);
hold on
c = tmul(c,u{3}',3);

style_est_results = zeros(numofletters,numofstyles);
letter_est_results = zeros(numofletters,numofstyles);
reconstruction_error = zeros(numofletters,numofstyles);
%for plot
plot_test_flg = 1;
estimated_letter_vector = zeros(numel(tSeqs),length(letter_bases));
estimated_letter_label = zeros(numel(tSeqs));
est_index = 1;
correct_letter_est = 0;
correct_style_est = 0;
for k = 1:1
    G2 = tmul(c,style_bases',2);%compute G2 for all style_bases;
    for li = 1:numofletters
        for si=1:numofstyles
            disp '------------------------------------------'
            fprintf('For letter %s \n', getLetterName(li,1.0));
            fprintf('For subject %s \n',getSubjectName(si));

            %lmean = mean()
            testSeq = tSeqs{li,k,si};
            cycleframenum = size(testSeq,2);
            %% step 1: embedding on unit circle of the test sequence and find mapping
            % learn mapping between the manifolds and the input space
            CF=learnmapping_grbf_regulize(testSeq');
            %% step : create vector b
            %B = CF*u{3}';
            b = reshape(CF,[],1);
            %% step: iterate for all persons and compute the letter vector
            %%corresponding to it and choose the peson with letter has min
            %%reconstruction error
            err_min = inf;

            if p_use_sampling
                G = tmul(G2,letter_bases,1);
                G_ = unfold(G,3);
                [Indx,Dist] = knnsearch(G_',b','K',length(style_bases)*length(letter_bases));
                %Dist is the distace between the vector b and every vector
                %in the matrix G.
                Dist(Indx) = Dist;
                D = reshape(Dist,length(letter_bases),length(style_bases));
                [m1,l_] = min(D);
                [m2,p] = min(m1);
                l = l_(p);
            else
                for i = 1:length(style_bases)
                    p_ = style_bases(i,:)';
                    %G2 = tmul(c,p_,2);%is to compute G2 each iteration is
                    %equivalent to take the ith slice of the 2nd dimension but
                    %with a small approx error, need investication 
                    %ToDO
                    % estimate corresponding letter vector
                    %l_ = unfold(G2(:,i,:),1)'\b;
                    l_ = squeeze(G2(:,i,:))'\b;
                    % estimate the reconstruction error
                    %b_bar = unfold(tmul(G2(:,i,:),l_,1),3);
                    b_bar = squeeze(G2(:,i,:))'*l_;
                    err = norm(b-b_bar);
                    if err<err_min
                        p = p_;
                        l = l_;
                        err_min = err;
                    end
                end
            end
            %p
            %l
            %b_bar = squeeze(tmul(tmul(c,l,1),p,2));
            final_error = err_min %norm(b-b_bar)
            reconstruction_error(li,si) = final_error;
            %% estimate the nearest latter and subject
            disp ' for plotting'
            estimated_letter_vector(est_index,:) = l';
            estimated_letter_label(est_index) = li;
            est_index = est_index+1;
            if(plot_test_flg)
                plotting(l(1:3),li,1,fig);
            end
            
            disp '**Results:'
            mltr = double(length(letter_bases))/numofletters;
            %use knnclassify instead of knnsearch, it does the
            %classificaton directly.
            if(p_use_wts)
                [ltr,D] = knnsearch(letter_bases,l','K',5,'Distance','seuclidean','Scale',ltr_wts);%?????????????
            else
                [ltr,D] = knnsearch(letter_bases,l','K',5);%?????????????
            end
            %ltr = ceil(ltr/mltr)%????????????
            ltr = mod(ltr-1,numofletters)+1;
            D
            ltr_unq = unique(ltr);%,'stable');
            ltr_rep = histc(ltr,ltr_unq)
            [dmy,idx] = max(ltr_rep)%ToDo if more than one max value then should choose based on sum of wights, currently chooose the lowest ranked one
            letter_est_results(li,si) = ltr_unq(idx);
            if ltr_unq(idx)==li
               correct_letter_est = correct_letter_est +1;
            end
            disp 'Nearest letter is:'
            getLetterName(ltr_unq(idx),1)
            %getLetterName(ltr,1.0)
            
            mstyle = double(length(style_bases))/numofstyles;
            [style,D] = knnsearch(style_bases,p');%,'K',10)
            style_est_results(li,si) =ceil(style/mstyle);%mod(style-1,numofstyles)+1;
            if ceil(style/mstyle)==si
                correct_style_est = correct_style_est+1;
            end
            disp 'Nearest Person is:'
            getSubjectName(ceil(style/mstyle))%mod(style-1,numofstyles)+1)
        end
    end
%     g = 1:numofstyles;
%     g = repmat(g,numofletters,1);
%     speaker_error_percent = nnz(style_est_results - g)/double(numel(g));
%     speaker_iden_rate = 1 - speaker_error_percent
    speaker_iden_rate = correct_style_est/numel(tSeqs)
%     g = 1:numofletters;
%     g = repmat(g',1,numofstyles);
%     letter_error_percent = nnz(letter_est_results - g)/double(numel(g))
    letter_recog_rate = correct_letter_est/numel(tSeqs)
end

%save('results_HoG-3-52-50','reconstruction_error','letter_est_results','style_est_results','letter_error_percent','speaker_error_percent');
%save('results_Gabor-set2_26-20-(1-2)','reconstruction_error','letter_est_results','style_est_results','letter_error_percent','speaker_error_percent');
save('results_52-10-(1-2)_NEW','reconstruction_error','letter_est_results','style_est_results','letter_recog_rate','speaker_iden_rate','estimated_letter_vector','estimated_letter_label');

diary off
return 

%%%%%%%%%%%%%%%%%%%%%%%%%%%
imagesize = 4800;
CONTENT_PRECISION = 50;
centerNum = 10;

NEARZERO = 0.000000000001;
MINSQUARE = 1;
for i = 1: numofstyles
    for j = 1: numofletters
        
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

            imagev = testSeq(:,k)';%FEIMAGE{i}{j}(k,:);
            
            
            
            %ec = mean(exp_mean,2);

%             Tpsi = tmul(Tp, style_mean(i,:)',1);
            Tpsc = tmul(Tp, pmean',2);
%             Tpsiec = tmul(Tpsi, ec',2);
            Tpscei = tmul(Tpsc, letter_bases(j,:)',1);            
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
    
