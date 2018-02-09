function [ Div ] = returnDiv( this, mode )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

    D = sparse(this.d0);
    rB = D*this.Mesh.X;
    rB = bsxfun(@rdivide,rB,this.Lp);

    DivX = bsxfun(@times,D,this.Ld.*rB(:,1));
    DivY = bsxfun(@times,D,this.Ld.*rB(:,2));
    DivZ = bsxfun(@times,D,this.Ld.*rB(:,3));

    Div1 = bsxfun(@times,DivX,this.e1(:,1)') + bsxfun(@times,DivY,this.e1(:,2)') + bsxfun(@times,DivZ,this.e1(:,3)');
    Div2 = bsxfun(@times,DivX,this.e2(:,1)') + bsxfun(@times,DivY,this.e2(:,2)') + bsxfun(@times,DivZ,this.e2(:,3)');
    if (nargin == 1 || strcmp(mode,'surface'))
        Div = sparse(-[Div1';Div2']);
    else
        n = bsxfun(@rdivide,this.K,sqrt(sum(this.K.^2,2)));
        Div3 = bsxfun(@times,DivX,n(:,1)') + bsxfun(@times,DivY,n(:,2)') + bsxfun(@times,DivZ,n(:,3)');
        Div = sparse(-[Div1';Div2';Div3']);
    end

end

