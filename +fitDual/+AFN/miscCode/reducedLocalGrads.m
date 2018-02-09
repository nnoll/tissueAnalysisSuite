function [ drX1, drY1, drX2, drY2 ] = reducedLocalGrads( q, p, bCells, d0, r1 ,r2 )
%REDUCEDLOCALGRADS Summary of this function goes here
%   Detailed explanation goes here

    Nb = size(d0,1);
    Nc = size(d0,2);
    Eind = (1:Nb)';

    rowX1 = zeros(4*Nb,1);
    colX1 = zeros(4*Nb,1);
    valX1 = zeros(4*Nb,1);
    
    rowX2 = zeros(4*Nb,1);
    colX2 = zeros(4*Nb,1);
    valX2 = zeros(4*Nb,1);

    rowY1 = zeros(4*Nb,1);
    colY1 = zeros(4*Nb,1);
    valY1 = zeros(4*Nb,1);

    rowY2 = zeros(4*Nb,1);
    colY2 = zeros(4*Nb,1);
    valY2 = zeros(4*Nb,1);
    
    n = 1;
    II = (n-1)*Nb + (1:Nb);
    
    rowX1(II) = Eind;
    colX1(II) = bCells(:,1);
    valX1(II) = p(bCells(:,1));
    
    rowX2(II) = Eind;
    colX2(II) = bCells(:,1);
    valX2(II) = p(bCells(:,1));
    
    rowY1(II) = Eind;
    colY1(II) = Nc + bCells(:,1);
    valY1(II) = p(bCells(:,1));
    
    rowY2(II) = Eind;
    colY2(II) = Nc + bCells(:,1);
    valY2(II) = p(bCells(:,1));
    
    n = 2;
    II = (n-1)*Nb + (1:Nb);
    
    rowX1(II) = Eind;
    colX1(II) = bCells(:,2);
    valX1(II) = -p(bCells(:,2));
    
    rowX2(II) = Eind;
    colX2(II) = bCells(:,2);
    valX2(II) = -p(bCells(:,2));
    
    rowY1(II) = Eind;
    colY1(II) = Nc + bCells(:,2);
    valY1(II) = -p(bCells(:,2));
    
    rowY2(II) = Eind;
    colY2(II) = Nc + bCells(:,2);
    valY2(II) = -p(bCells(:,2));
    
    n = 3;
    II = (n-1)*Nb + (1:Nb);

    rowX1(II) = Eind;
    colX1(II) = 2*Nc + bCells(:,1);
    valX1(II) = q(bCells(:,1),1) - r1(:,1);
    
    rowX2(II) = Eind;
    colX2(II) = 2*Nc + bCells(:,1);
    valX2(II) = q(bCells(:,1),1) - r2(:,1);
    
    rowY1(II) = Eind;
    colY1(II) = 2*Nc + bCells(:,1);
    valY1(II) = q(bCells(:,1),2) - r1(:,2);
    
    rowY2(II) = Eind;
    colY2(II) = 2*Nc + bCells(:,1);
    valY2(II) = q(bCells(:,1),2) - r2(:,2);

    n = 4;
    II = (n-1)*Nb + (1:Nb);

    rowX1(II) = Eind;
    colX1(II) = 2*Nc + bCells(:,2);
    valX1(II) = r1(:,1) - q(bCells(:,2),1);
    
    rowX2(II) = Eind;
    colX2(II) = 2*Nc + bCells(:,2);
    valX2(II) = r2(:,1) - q(bCells(:,2),1); 
    
    rowY1(II) = Eind;
    colY1(II) = 2*Nc + bCells(:,2);
    valY1(II) = r1(:,2) - q(bCells(:,2),2);
    
    rowY2(II) = Eind;
    colY2(II) = 2*Nc + bCells(:,2);
    valY2(II) = r2(:,2) - q(bCells(:,2),2);
   
    drX1 = sparse(rowX1,colX1,valX1,Nb,3*Nc);
    drX2 = sparse(rowX2,colX2,valX2,Nb,3*Nc);
    drY1 = sparse(rowY1,colY1,valY1,Nb,3*Nc);
    drY2 = sparse(rowY2,colY2,valY2,Nb,3*Nc);

end

%     drX1 = zeros(Nb,3*Nc); % Sparse gradient of unnormed vector.
%     drY1 = zeros(Nb,3*Nc); % Sparse gradient of unnormed vector.
% 
%     drX2 = zeros(Nb,3*Nc); % Sparse gradient of unnormed vector.
%     drY2 = zeros(Nb,3*Nc); % Sparse gradient of unnormed vector.

    % X gradient component
%     drX1(Eind + Nb*(bCells(:,1)-1)) = p(bCells(:,1)); 
%     drX1(Eind + Nb*(bCells(:,2)-1)) = -p(bCells(:,2)); 
%     drX2(Eind + Nb*(bCells(:,1)-1)) = p(bCells(:,1)); 
%     drX2(Eind + Nb*(bCells(:,2)-1)) = -p(bCells(:,2));
% 
%     % Y gradient component
%     drY1(Eind + Nb*(Nc + bCells(:,1)-1)) = p(bCells(:,1)); 
%     drY1(Eind + Nb*(Nc + bCells(:,2)-1)) = -p(bCells(:,2)); 
%     drY2(Eind + Nb*(Nc + bCells(:,1)-1)) = p(bCells(:,1)); 
%     drY2(Eind + Nb*(Nc + bCells(:,2)-1)) = -p(bCells(:,2));
% 

%     % P gradient component
%     drX1(Eind + Nb*(2*Nc + bCells(:,1)-1)) =  q(bCells(:,1),1) - r1(:,1); 
%     drX1(Eind + Nb*(2*Nc + bCells(:,2)-1)) = r1(:,1) - q(bCells(:,2),1); 
%     drY1(Eind + Nb*(2*Nc + bCells(:,1)-1)) = q(bCells(:,1),2) - r1(:,2); 
%     drY1(Eind + Nb*(2*Nc + bCells(:,2)-1)) = r1(:,2) - q(bCells(:,2),2);
% 
%     drX2(Eind + Nb*(2*Nc + bCells(:,1)-1)) = q(bCells(:,1),1) - r2(:,1); 
%     drX2(Eind + Nb*(2*Nc + bCells(:,2)-1)) = r2(:,1) - q(bCells(:,2),1); 
%     drY2(Eind + Nb*(2*Nc + bCells(:,1)-1)) = q(bCells(:,1),2) - r2(:,2); 
%     drY2(Eind + Nb*(2*Nc + bCells(:,2)-1)) = r2(:,2) - q(bCells(:,2),2);
% 
%     drX1 = sparse(drX1);
%     drY1 = sparse(drY1);
%     drX2 = sparse(drX2);
%     drY2 = sparse(drY2);
