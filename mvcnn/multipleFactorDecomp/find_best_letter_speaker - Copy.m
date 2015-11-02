function [l,p,err]=find_best_letter_speaker(b,c,init_value)
        threshold = 1E-4;
        b_bar_old = zeros(size(b),1);
        
        err = 1.0;err_old = 2.0;
        count = 1;
        l = init_value;
        while (err> threshold) && (err_old>=err)
            % estimate corresponding person style vector
            G1 = tmul(c,l,1);
            p = unfold(G1,2)'\b;
%             %redirect p vector in the column space of style_bases
% doesn't have meaning because the bases are orthonotmal so any vector will
% by default belong to the there space. this will be meaningful if the
% bases are not independent.
%             x = style_bases'\p;
%             p = style_bases'*x;
            % estimate corresponding letter vector
            G2 = tmul(c,p,2);
            l = unfold(G2,1)'\b;
%             %redirect p vector in the column space of style_bases
%             x = letter_bases'\l;
%             l = letter_bases'*x;

            b_bar = unfold(tmul(G2,l,1),3);

            err = norm(b_bar-b_bar_old);
            %error = norm(b-b_bar)
            b_bar_old = b_bar;
            err_old = err;
            count = count+1
            if count>20
                break;
            end
        end
end