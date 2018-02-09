function [ Grad ] = returnGrad( this, mode )

    if (nargin == 1 || strcmp(mode,'2D'))
        D = this.d0;
        rB = D*this.Mesh.X;
        rB = bsxfun(@rdivide,rB,this.Lp.^2);

        GradX = bsxfun(@times,D,rB(:,1));
        GradY = bsxfun(@times,D,rB(:,2));
        GradZ = bsxfun(@times,D,rB(:,3));

        Grad1 = bsxfun(@times,GradX,this.e1(:,1)') + bsxfun(@times,GradY,this.e1(:,2)') + bsxfun(@times,GradZ,this.e1(:,3)');
        Grad2 = bsxfun(@times,GradX,this.e2(:,1)') + bsxfun(@times,GradY,this.e2(:,2)') + bsxfun(@times,GradZ,this.e2(:,3)');
        Grad = [Grad1,Grad2];
    else
        D = this.d0;
        rB = D*this.Mesh.X;
        rB = bsxfun(@rdivide,rB,this.Lp.^2);

        GradX = bsxfun(@times,D,rB(:,1));
        GradY = bsxfun(@times,D,rB(:,2));
        GradZ = bsxfun(@times,D,rB(:,3));

        Grad = [GradX,GradY,GradZ];
    end

end

