function [ r ] = computePrimalVerts( this, rV )
    % COMPUTE PRIMAL VERTS. Given the stored dual triangulation, this
    % function will compute and return the position of the vertices in the
    % underlying primal cell decomposition.
    
    if (var(this.p) > 0)
        % Initialize rotation matrix.
        Rot = [0,-1;1,0];
        tri = this.pTri;

        % Store pressure, theta, and q positions.
        p1 = this.p(tri(:,1));
        p2 = this.p(tri(:,2));
        p3 = this.p(tri(:,3));

        z1 = this.theta(tri(:,1));
        z2 = this.theta(tri(:,2));
        z3 = this.theta(tri(:,3));

        q1 = this.q(tri(:,1),:);
        q2 = this.q(tri(:,2),:);
        q3 = this.q(tri(:,3),:);

        % Store triangular edges and its area.
        t1 = bsxfun(@times,p3,q3) - bsxfun(@times,p2,q2);
        t2 = bsxfun(@times,p1,q1) - bsxfun(@times,p3,q3);
        t3 = bsxfun(@times,p2,q2) - bsxfun(@times,p1,q1);

        S0 = .5*abs(t1(:,1).*t2(:,2) - t1(:,2).*t2(:,1));

        r0 =  (bsxfun(@times,(p1.*(sum(q1.^2,2)) + z1),t1) ... 
                + bsxfun(@times,(p2.*(sum(q2.^2,2)) + z2),t2) ...
                + bsxfun(@times,(p3.*(sum(q3.^2,2)) + z3),t3))*Rot';
        r0 = bsxfun(@rdivide,r0,4*S0);

        % Compute the deviation off the zeroth order term.
        eta = bsxfun(@times,p1,t1) + bsxfun(@times,p2,t2) + bsxfun(@times,p3,t3);
        etaMag = sqrt(sum(eta.^2,2));
        etaStar = eta*Rot';
        etaStar = bsxfun(@rdivide,etaStar,etaMag);

        delta = ((2*S0./etaMag) - dot(etaStar,r0,2)) - sqrt( (((2*S0./etaMag) - dot(etaStar,r0,2))).^2 - sum(r0.^2,2));
        delta2 = ((2*S0./etaMag) - dot(etaStar,r0,2)) + sqrt( (((2*S0./etaMag) - dot(etaStar,r0,2))).^2 - sum(r0.^2,2));

        delta(isnan(delta)) = 0;
        delta = real(delta);

        delta2(isnan(delta2)) = 0;
        delta2 = real(delta2);
        
        if (nargin == 2)
            etaStar(isnan(etaStar)) = 0;
            deltaR1 = bsxfun(@times,delta,etaStar);
            deltaR2 = bsxfun(@times,delta2,etaStar);

            r1 = r0 + deltaR1;
            r2 = r0 + deltaR2;

            D1 = sum((r1-rV).^2,2);
            D2 = sum((r2-rV).^2,2);
            [~,ind] = min([D1,D2],[],2);

            r = r1;
            r(ind==2,:) = r2(ind==2,:);
        else
            ind = sum(this.tri == this.pTri,2);
            ind = (ind == 1);
            delta(ind) = delta2(ind);

            etaStar(isnan(etaStar)) = 0;
            deltaR = bsxfun(@times,delta,etaStar);
            r = r0 + deltaR;
        end
        
    else
        
        Rot = [0,-1;1,0];
        tri = this.tri;

        % Store theta, and q positions.
        z1 = this.theta(tri(:,1));
        z2 = this.theta(tri(:,2));
        z3 = this.theta(tri(:,3));

        q1 = this.q(tri(:,1),:);
        q2 = this.q(tri(:,2),:);
        q3 = this.q(tri(:,3),:);

        % Store triangular edges and its area.
        t1 = q3 - q2;
        t2 = q1 - q3;
        t3 = q2 - q1;

        S0 = .5*abs(t1(:,1).*t2(:,2) - t1(:,2).*t2(:,1));

        r =  (bsxfun(@times,(sum(q1.^2,2) + z1),t1) ... 
            + bsxfun(@times,(sum(q2.^2,2) + z2),t2) ...
            + bsxfun(@times,(sum(q3.^2,2) + z3),t3))*Rot';
        r = bsxfun(@rdivide,r,4*S0);
        
    end
    
end

