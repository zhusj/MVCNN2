function [l,p,err]=find_best_letter_speaker_ltr_sampling(b,Core,letter_bases,style_bases)
    new_letter_bases = letter_bases;%every row represent one bases
    for ewer = 1:5
        G = tmul(G2,new_letter_bases',1); %verified
        G_ = reshape(G,size(new_letter_bases,1)*size(style_bases,1),[]);%unfold(G,3);
        Dist = pdist2(G_,b'); %[Indx,Dist] = knnsearch(G_',b','K',size(new_letter_bases,1)*size(style_bases,1));%verified
        %Dist is the distace between the vector b and every vector
        %in the matrix G.
        %Dist(Indx) = Dist;%verified
        D = reshape(Dist,size(new_letter_bases,1),size(style_bases,1));%verified
        [mL,l_] = min(D);%verified
        [mP,p_] = min(mL);%veririfed
        p = style_bases(p_,:)';%verified
        best_style = p_;
        l = new_letter_bases(l_(p_),:)';%verified
        break;
        %
        b_bar = squeeze(G2(:,p_,:))'*l;
        err = norm(b-b_bar)%err_min 
        if(err<.8)
            break;
        end
        %resample
        d = D(:,p_);%min(D,[],2);
        [d,indx] = sort(d);
        new_letter_bases = resample_new_bases(new_letter_bases,d,indx);%verified
        if size(new_letter_bases,1)==1
            l = new_letter_bases';
            break;
        end
    end
end