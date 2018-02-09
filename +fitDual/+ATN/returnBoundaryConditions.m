function [ Aeq, beq ] = returnBoundaryConditions( x0, bulkCells, extCells, extBonds, rBext )
    % RETURN BOUNDARY CONDITIONS 

    Aeq = zeros(size(extBonds,1),length(x0));
    extCells = (length(bulkCells)+1):(length(bulkCells)+length(extCells)); 

    C = length(x0)/3;

    for b = 1:size(extBonds,1)
       Aeq(b,extCells(extBonds(b,1))) = rBext(b,1);
       Aeq(b,extCells(extBonds(b,2))) = -rBext(b,1);
       Aeq(b,C+extCells(extBonds(b,1))) = rBext(b,2);
       Aeq(b,C+extCells(extBonds(b,2))) = -rBext(b,2);
    end

    beq = zeros(size(extBonds,1),1);
       
end

