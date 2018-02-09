function [ L, bndryCells ] = returnSmoothingOperator( this )
%RETURNSMOOTHINGOPERATOR Summary of this function goes here
%   Detailed explanation goes here

    L = -this.d0'*this.d0;
    bndryEdges = sum(abs(this.d1),1) == 1;
    bndryCells = sum(abs(this.d0(bndryEdges,:)),1) >=1;
    L(bndryCells,:) = [];
    L(:,bndryCells) = [];

    L = L - diag(sum(L,1));
    L = .5*L;

end

