% tS = 25;

clear Tg c
for t = timePts
    t
    [ T, P ] = MI.invertMech(tStruct(t),1,1);
    [ i1 ] = generate.bondMap( tStruct(t) );
    Tg{t} = zeros(size(T));
    for b = 1:length(i1{1})
        if (i1{1}(b) > 0)
            Tg{t}(b) = tStruct(t).Bdat(i1{1}(b)).tension;
        else
            Tg{t}(b) = 0;
        end
    end

    % scatter(Tg{t}(Tg{t}>0),T(Tg{t}>0))
    c(t) = corr(Tg{t}(Tg{t}>0),T(Tg{t}>0));
end