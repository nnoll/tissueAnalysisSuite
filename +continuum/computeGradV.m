function [ vGrad, xG, yG ] = computeGradV( x, y, v )
    % COMPUTE GRAD V

    LY = (length(y)-1);
    LX = (length(x)-1);
    vGrad = cell(length(v),1);
    
    xG = .5*(x(1:(end-1)) + x(2:end));
    yG = .5*(y(1:(end-1)) + y(2:end));

    for t = 1:length(v)
        vGrad{t} = zeros(2,2,LY,LX);
        for ii = 1:LY
            for jj = 1:LX
                v1 = v{t}(:,ii,jj);
                v2 = v{t}(:,ii+1,jj);
                v3 = v{t}(:,ii,jj+1);
                v4 = v{t}(:,ii+1,jj+1);

                b = [v1(1);v2(1);v3(1);v4(1);v1(2);v2(2);v3(2);v4(2)];
                dx = [x(jj)-xG(jj);x(jj)-xG(jj);x(jj+1)-xG(jj);x(jj+1)-xG(jj)];
                dy = [y(ii)-yG(ii);y(ii+1)-yG(ii);y(ii)-yG(ii);y(ii+1)-yG(ii)];

                c1 = [ones(4,1);zeros(4,1)];
                c2 = [zeros(4,1);ones(4,1)];
                c3 = [dx;zeros(4,1)];
                c4 = [dy;zeros(4,1)];
                c5 = [zeros(4,1);dx];
                c6 = [zeros(4,1);dy];

                L = [c1,c2,c3,c4,c5,c6];
                params = L \ b;
                vGrad{t}(1,1,ii,jj) = params(3);
                vGrad{t}(1,2,ii,jj) = params(4);
                vGrad{t}(2,1,ii,jj) = params(5);
                vGrad{t}(2,2,ii,jj) = params(6);
            end 
        end
    end

%     [xG,yG] = meshgrid(xG,yG);

end

