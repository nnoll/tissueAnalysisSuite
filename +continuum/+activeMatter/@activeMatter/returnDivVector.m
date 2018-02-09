function [ divV ] = returnDivVector( this )
    % RETURN DIV VECTOR 

    % Use Green's first identity
    gradS = this.returnGradScalar();
    divV = -gradS';
    
end

