function [l,p,err]=find_best_letter_speaker(b,Core,letter_bases,style_bases)
    err_min = inf;
    if p_use_sampling
        G = tmul(Core,letter_bases,1);
        G_ = unfold(G,3);
        [Indx,Dist] = knnsearch(G_',b','K',length(style_bases)*length(letter_bases));
        %Dist is the distace between the vector b and every vector
        %in the matrix G.
        Dist(Indx) = Dist;
        D = reshape(Dist,length(letter_bases),length(style_bases));
        [m1,l_] = min(D);
        [m2,p_] = min(m1);
        p = style_bases(p_,:)';
        l = letter_bases(l_(p_),:)';
    else
        for i = 1:length(letter_bases)
            l_ = mean(letter_bases(i,:))';

            %[l_,p_,err]=find_best_letter_speaker(b,c,l_);
            %G2 = tmul(c,p_,2);%is to compute G2 each iteration is
            %equivalent to take the ith slice of the 2nd dimension but
            %with a small approx error, need investication 
            %ToDO
            % estimate corresponding letter vector
            %l_ = unfold(G2(:,i,:),1)'\b;


            p_ = squeeze(Core(i,:,:))'\b;

            %b_bar = tmul(c,l_,1);
            %b_bar = squeeze(b_bar)'*p_;
            % estimate the reconstruction error
            %b_bar = unfold(tmul(G2(:,i,:),l_,1),3);
            b_bar = squeeze(Core(i,:,:))'*p_;
            err = norm(b-b_bar);

            if err<err_min
                p = p_;
                l = l_;
                best_letter = i;
                err_min = err;
            end
        end
    end
end