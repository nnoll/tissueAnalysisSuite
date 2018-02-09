function [ Phi ] = symLR( Phi )
    % SYM LR 

   for t = 1:length(Phi)
        dim = size(Phi{t});
        if (length(dim) == 3)
            if (mod(size(Phi{t},2),2) == 0)
                L = size(Phi{t},2)/2;
                Lend = size(Phi{t},2);
                delta = 1;
            else
                L = (size(Phi{t},2)+1)/2;
                Lend = size(Phi{t},2);
                delta = 0;
            end
            Phi{t}(1,1:L,:) = .5*(Phi{t}(1,1:L,:) + Phi{t}(1,Lend:-1:(L+delta),:));
            Phi{t}(1,Lend:-1:(L+delta),:) = Phi{t}(1,1:L,:);
            Phi{t}(2,1:L,:) = .5*(Phi{t}(2,1:L,:) - Phi{t}(2,Lend:-1:(L+delta),:));
            Phi{t}(2,Lend:-1:(L+delta),:) = -Phi{t}(2,1:L,:);
        else

        end
   end
    
end

