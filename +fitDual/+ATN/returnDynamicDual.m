function [ PN, ERes ] = returnDynamicDual( Struct )
    %RETURN DYNAMIC DUAL 
    
    for t = 1:length(Struct)
        t
        if (t == 1)
            [ q, theta, tri, cells, ERes(t) ] = fitATN.returnDual( Struct(t) );
            PN{1} = pressure.net(q,theta,ones(size(theta)),tri,cells);
        else
            [ q, theta, tri, cells, ERes(t) ] = fitATN.returnDual( Struct(t) );
            PN{t} = pressure.net(q,theta,ones(size(theta)),tri,cells);
        end
    end

end

