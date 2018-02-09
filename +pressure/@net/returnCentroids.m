function [ rC ] = returnCentroids( this )

        [ faces ] = this.computeFaces( this.computePrimalVerts() );
        rV = this.computePrimalVerts();

        rC = zeros(size(faces,1),2);
        for c = 1:size(faces,1)
            verts = faces(c,~isnan(faces(c,:)));
            rC(c,:) = mean(rV(verts,:),1);
        end
end

