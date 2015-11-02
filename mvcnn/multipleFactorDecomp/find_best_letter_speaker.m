function [l,p,err]=find_best_letter_speaker(b,c,init_value)
        threshold = 1E-6;
        b_bar = zeros(size(b),1);
        sz = size(c);
        err = 1.0;err_old = 2.0;
        count = 1;
        l = init_value;
        c_1 = unfold(c,1);
        c_2 = unfold(c,2);
        while 1%(abs(err-err_old)> threshold) && (err_old>err)
            % estimate corresponding person style vector
            %tmul(c,l,1);
            G1 = reshape(l'*c_1,sz(2),[]);
            p = G1'\b;%unfold(G1,2)'\b;
            
            %error = norm(b-b_bar)
            b_bar_old = b_bar;
            %err_old = err;
            
            b_bar = (p'*G1)';%unfold(tmul(G2,l,1),3);
            err = norm(b_bar-b_bar_old)
            
            if err< threshold%(abs(err-err_old)< threshold)% || (err>err_old)
                break;
            end
%             %redirect p vector in the column space of style_bases
% doesn't have meaning because the bases are orthonotmal so any vector will
% by default belong to the there space. this will be meaningful if the
% bases are not independent.
%             x = style_bases'\p;
%             p = style_bases'*x;
            % estimate corresponding letter vector
            G2 = reshape(p'*c_2,sz(1),[]);
            %G2 = tmul(c,p,2);
            l = G2'\b;
%             %redirect p vector in the column space of style_bases
%             x = letter_bases'\l;
%             l = letter_bases'*x;

            
            
            count = count+1
            if count>20
                break;
            end
        end
        err = norm(b_bar-b);
 end