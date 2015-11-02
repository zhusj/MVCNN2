function [l,p,err]=get_best_sampling(b,Core,D1_bases,D2_bases)
    thresh = 0.001;
    letter_samples = D1_bases;%every row represent one bases
    style_samples = D2_bases;
    for ewer = 1:20
        G = tmul(tmul(Core,letter_samples',1),style_samples',2); %verified
        G_ = reshape(G,size(letter_samples,1)*size(style_samples,1),[]);%unfold(G,3);
        Dist = dist2(G_,b'); %[Indx,Dist] = knnsearch(G_',b','K',size(new_letter_bases,1)*size(style_bases,1));%verified
        %Dist is the distace between the vector b and every vector
        %in the matrix G.
        %Dist(Indx) = Dist;%verified
        D = reshape(Dist,size(letter_samples,1),size(style_samples,1));%verified
        
        [mL,p_list] = min(D,[],2);%verified
        [~,l_i] = min(mL);%get the over all minimum value in D, l_ will be the letter index
        p_i = p_list(l_i);%l_, p_ will be the index of the over all minimum value in D
        %p = style_bases(p_,:)';%verified
        %best_style = p_;
        %l = letter_samples(l_(p_),:)';%verified
        %break;
        %
        
        l = letter_samples(l_i,:)';
        p = style_samples(p_i,:)';
        %b_bar = %queeze(G2(:,p_,:))'*l;
        err = D(l_i,p_i);%norm(b-b_bar)%err_min 
        if(err<thresh)
            break;
        end
        %if threshold has not been met then resample new letter and speaker
        %samples based on weights represented by mimimum of all rows and
        %columns.
        
        letter_samples = generate_new_samples(letter_samples(index,:),mL);%verified
        
        mP = min(D);%veririfed
        [d,indx] = sort(mP);
        style_samples = generate_new_samples(style_samples,d,indx);%verified
        
        
        if size(letter_samples,1)==1
            l = letter_samples';
            break;
        end
    end
    ewer
end