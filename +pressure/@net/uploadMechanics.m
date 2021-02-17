function [ Struct, found_bond ] = uploadMechanics( this, Struct )
    % uploadMechanics Store details of the mechanics of the pressure net in
    % Struct
    %
    % Parameters
    % ----------
    % this : PN{t} instance -- a pressure net at a given time 
    %        Do not supply this argument, since it is implied from
    %        uploadMechanics being a method of net
    % Struct : the data structure
    % 
    % Returns
    % -------
    % Struct : data structure
    %   contains the following fields:
    % vdat: 
    %     nverts : neighbor list for vertices
    %     ncells : 
    %     vertxcoord : column of data in which vertex lives
    %     vertycoord : row of data in which each vertex lives
    % cdat : cell data    
    %     ncells : indices of neighboring cells
    %     nverts : vertices that define the cell
    %     centroid.coord(1) x position of cell centroid
    %     centroid.coord(2) y position of cell centroid
    % bdat : 
    %     nverts
    %     ncells
    %     pix : linear indices of the pixels associated with that bond

    % found_bonds : list
    %   list of indices for which we have found a matching bond

    [ T ] = this.returnTension();
    [ p ] = this.p;
    
    found_bond = zeros(length(Struct.Bdat),1);
    % iterate through bonds
    for b = 1:length(Struct.Bdat)
        bCells = Struct.Bdat(b).cells;
        % find the cells that the bond corresponds to 
        if (sum(ismember(bCells,this.cellLabels)) == 2) % Both cells involved in fit.
            % Find indices of the cells where this bond row 
            cellInd1 = find(ismember(this.cellLabels,bCells(1)));
            cellInd2 = find(ismember(this.cellLabels,bCells(2)));
            % d0 is an e x c matrix of exterior derivatives with +1 and -1s
            % at the endpts of each bond
            % d1 is a v x e matrix of exterior derivatives. Upstream is +1,
            % downstream is -1 when moving counterclockwise around a
            % tension plaquette.
            % d0 and d1 are matrices that take derivatives 
            bInd = find(this.d0(:, cellInd1)~=0 & this.d0(:, cellInd2)~=0);
            
            % Enforce the condition that if we are missing tension, do not
            % include chem, which is the integrated intensity per edge for
            % that fluorescence channel (that 'chem')
            if (isfield(Struct.Bdat(b), 'chem'))
                if (~isempty(bInd) && ~isempty(Struct.Bdat(b).chem))
                    Struct.Bdat(b).tension = T(bInd);
                    found_bond(b) = 1;
                end
            else
                if (~isempty(bInd))
                    Struct.Bdat(b).tension = T(bInd);
                    found_bond(b) = 1;
                end
            end
        end
    end
    
    found_bond = find(found_bond);
    
    % insert pressure into Struct.Cdat().pressure
    for c = 1:length(p)
        cellLabel = this.cellLabels(c);
        Struct.Cdat(cellLabel).pressure = p(c);
    end
    
end

