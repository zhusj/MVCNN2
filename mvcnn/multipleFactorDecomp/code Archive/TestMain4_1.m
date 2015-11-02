function TestMain4_1
% FILENAME: tensor_facial_exp_eval1.m
% This experimnet for speaker identification
%train over 26x2x10, test for 26x1x10.

clear
pack
diary('run_log_3.txt');
%load('TEST_DATA.mat')% 'c','u','s','rdim'
disp 'train over reptns 1,2 and test over reptn 3'
subjects = {'*'}
letters = {'*'}
repeats = {'3'}
[tSeqs,tSizes,tLetters,tSubjects] = LoadDataAVLetterFiles(letters,subjects,repeats);


% % load test data
%load('AVLetter_Train_SVD_26-10-2_reps-mean');
load('AVLetter_Train_SVD_26-10-(1-2)_merge-letters-reps');
%load('AVLetter_Train_SVD-merge_letters_reps');
%load('AVLetter_Train_SVD_merge-reps-subjs');
%c core tensor
%u{1} letter bases
%u{2} person style bases
threshold = 1E-10;
numofletters = size(tSeqs,1);%size(c,1);
numofstyle = size(tSeqs,3);%size(c,2);
pixelDIM = size(c,3);
letter_bases = u{1};
style_bases = u{2};
% %compute the similarity matrices for the bases
% let_bs = letter_bases/abs(max(max(letter_bases)));
% letter_similarity = let_bs'*let_bs;
% stl_bs = style_bases/abs(max(max(style_bases)));
% style_similarity = stl_bs'*stl_bs;
% save('similarity_mats','letter_similarity','style_similarity');
%Tp =tmul(c, u{3}',3);
style_mean = mean(style_bases,1)';% estimate the initial of person style vectors from the trained style basis vectors 
letter_mean = mean(letter_bases,1)';
c = tmul(c,u{3}',3);

style_est_results = zeros(numofletters,numofstyle);
letter_est_results = zeros(numofletters,numofstyle);
reconstruction_error = zeros(numofletters,numofstyle);

for k = 1:1
    %G1 = tmul(c,letter_bases',1);%compute G2 for all style_bases;
    for li = 1:numofletters
        for si=1:numofstyle
            disp '------------------------------------------'
            fprintf('For letter %s \n', getLetterName(li,numofletters/26.0));
            fprintf('For subject %s \n',getSubjectName(si));

            %lmean = mean()
            testSeq = featureExtraction(tSeqs{li,k,si});
            cycleframenum = size(testSeq,2);
            %% step 1: embedding on unit circle of the test sequence and find mapping
            disp('learn mapping between the manifolds and the input space');
            frames=cycleframenum;
            t=[1:frames]'*2*si/frames;
            P=[cos(t) sin(t)];
            % learn mapping between the manifolds and the input space
            CF=learnmapping_grbf(testSeq',P);
            %% step : create vector b
            %B = CF*u{3}';
            b = reshape(CF,[],1);
            %% step: iterate for all persons and compute the letter vector
            %%corresponding to it and choose the peson with letter has min
            %%reconstruction error
            err_min = inf;

            
            for i = 1:length(letter_bases)
                l_ = letter_bases(i,:)';
                G1 = tmul(c,l_,1);
                %G2 = tmul(c,p_,2);%is to compute G2 each iteration is
                %equivalent to take the ith slice of the 2nd dimension but
                %with a small approx error, need investication 
                %ToDO
                % estimate corresponding letter vector
                p_ = unfold(G1,2)'\b;
                % estimate the reconstruction error
                b_bar = unfold(tmul(G1,p_,2),3);
                err = norm(b-b_bar);
                if err<err_min
                    p = p_;
                    l = l_;
                    err_min = err;
                end
            end
            %p
            %l
            %b_bar = squeeze(tmul(tmul(c,l,1),p,2));
            final_error = err_min %norm(b-b_bar)
            reconstruction_error(li,si) = final_error;
            %% estimate the nearest latter and subject
            disp '**Results:'
            [ltr,D] = knnsearch(letter_bases,l','K',5);%?????????????
            ltr = ceil(ltr/2.0)
            D
            ltr_unq = unique(ltr);%,'stable');
            ltr_rep = histc(ltr,ltr_unq)
            [dmy,idx] = max(ltr_rep)%ToDo if more than one max value then should choose based on sum of wights, currently chooose the lowest ranked one
            letter_est_results(li,si) = ltr_unq(idx);
            disp 'Nearest letter is:'
            getLetterName(ltr,1.0)
            
            [style,D] = knnsearch(style_bases,p')%,'K',10)
            style_est_results(li,si) = style(1);
            disp 'Nearest Person is:'
            getSubjectName(style)
        end
    end
    g = 1:numofstyle;
    g = repmat(g,numofletters,1);
    speaker_error_percent = nnz(style_est_results - g)/double(numel(g))
    g = 1:numofletters;
    g = repmat(g',1,numofstyle);
    letter_error_percent = nnz(letter_est_results - g)/double(numel(g))
end
save('results_3','reconstruction_error','letter_est_results','style_est_results','letter_error_percent','speaker_error_percent');
diary off
return 

%%%%%%%%%%%%%%%%%%%%%%%%%%%
imagesize = 4800;
CONTENT_PRECISION = 50;
centerNum = 10;

NEARZERO = 0.000000000001;
MINSQUARE = 1;
for i = 1: numofstyle
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
    
