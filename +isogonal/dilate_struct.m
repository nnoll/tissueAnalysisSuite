function [ cStruct, Theta, three_cells ] = dilate_struct( Struct, a, L )
%DILATE_STRUCT - Using a generated reference lattice - this code will
%`dilate' each cell by a random number to create a fluctuating lattice.
%Will return the new vertex positions as well as the random numbers used. 

%%Inputs
%1. Struct - Generated Struct
%2. Alpha - Dilation variance.
%3. Beta - Randomness toggle

%%Outputs
%1. cStruct - Fluctuating Struct
%2. Theta - Parameters of fluctuation.

cStruct(1) = Struct;

%% Generate the forward matrix

%Find the verts that are only connected to internal three cells.
C = length(Struct.Cdat);
cells = 1:C;
three_fold = [Struct.Cdat.all_threefold];
three_cells = cells(find(three_fold));
% in_view = [Struct.Cdat.inview];
% view_cells = cells(find(in_view));
% three_cells = three_cells(ismember(three_cells,view_cells));

three_verts = [];
for tc = three_cells
    three_verts = horzcat(three_verts,Struct.Cdat(tc).nverts);
end

[uvert,ncounts] = count_unique(three_verts);
three_verts = uvert;

F = zeros(2*length(three_verts),length(three_cells));
for ii=1:length(three_cells)
    cell = three_cells(ii);
    cverts = Struct.Cdat(cell).orderednverts(1:end-1);
    cangle = zeros(2,length(cverts));
    cind = zeros(1,length(cverts));
    n_ext = zeros(2,length(cverts));
    for jj=1:length(cverts)
        cv = cverts(jj);
        %Get external vert
        nverts = Struct.Vdat(cv).nverts;
        ev = nverts(~ismember(nverts,cverts));
        n_ext(:,jj) = [Struct.Vdat(ev).vertxcoord - Struct.Vdat(cv).vertxcoord; Struct.Vdat(ev).vertycoord-Struct.Vdat(cv).vertycoord];
        n_ext(:,jj) = n_ext(:,jj)/norm(n_ext(:,jj));
        %Find index of vert
        if (ismember(cv,three_verts))
            cind(jj) = find(three_verts==cv,1);
        end
        %Get angles
        if (jj == 1)
            vm = [Struct.Vdat(cverts(length(cverts))).vertxcoord-Struct.Vdat(cv).vertxcoord;Struct.Vdat(cverts(length(cverts))).vertycoord-Struct.Vdat(cv).vertycoord];
            vp = [Struct.Vdat(cverts(jj+1)).vertxcoord-Struct.Vdat(cv).vertxcoord;Struct.Vdat(cverts(jj+1)).vertycoord-Struct.Vdat(cv).vertycoord];
        elseif (jj == length(cverts))
            vm = [Struct.Vdat(cverts(jj-1)).vertxcoord-Struct.Vdat(cv).vertxcoord;Struct.Vdat(cverts(jj-1)).vertycoord-Struct.Vdat(cv).vertycoord];
            vp = [Struct.Vdat(cverts(1)).vertxcoord-Struct.Vdat(cv).vertxcoord;Struct.Vdat(cverts(1)).vertycoord-Struct.Vdat(cv).vertycoord];
        else
            vm = [Struct.Vdat(cverts(jj-1)).vertxcoord-Struct.Vdat(cv).vertxcoord;Struct.Vdat(cverts(jj-1)).vertycoord-Struct.Vdat(cv).vertycoord];
            vp = [Struct.Vdat(cverts(jj+1)).vertxcoord-Struct.Vdat(cv).vertxcoord;Struct.Vdat(cverts(jj+1)).vertycoord-Struct.Vdat(cv).vertycoord];
        end
        
        vm = vm/norm(vm); vp = vp/norm(vp);
        angle(2) = atan2(norm(cross(vertcat(vm,0),vertcat(n_ext(:,jj),0))),dot(vertcat(vm,0),vertcat(n_ext(:,jj),0)));
        angle(1) = atan2(norm(cross(vertcat(vp,0),vertcat(n_ext(:,jj),0))),dot(vertcat(vp,0),vertcat(n_ext(:,jj),0)));
        cangle(:,jj) = sin(angle);
    end
       
    % Use geometric factors to build forward matrix.
    for jj=1:length(cverts)
        if (ismember(cverts(jj),three_verts))
            clear shift_cangle
            shift_cangle(1,:) = circshift(cangle(1,:)',size(cangle,2)-jj+1)';
            shift_cangle(2,:) = circshift(cangle(2,:)',size(cangle,2)-jj)';
            shift_cangle(:,size(shift_cangle,2)) = [];
            clear f
            f = 0;
            for kk=1:size(shift_cangle,2)
                clear g
                g = 1;
                for ll=1:kk
                    g = g * (shift_cangle(1,ll)/shift_cangle(2,ll));
                end
                f = f + g;
            end
            
            f = 1 / (1+f);
            F(2*cind(jj)-1,ii) = n_ext(1,jj)*f;
            F(2*cind(jj),ii) = n_ext(2,jj)*f;
        end
    end
end

%% Generate dilation parameters. Use to calculate perturbations.
Theta = a*rand(length(three_cells),1);
% X0 = 425;
% for ii=1:length(three_cells)
%     c = three_cells(ii);
%     Theta(ii) = a*(1-2*exp(-((Struct(1).Cdat(c).centroid.coord(1)-X0)/(3*sqrt(3)*L)).^2));    
% end
dr = F*Theta;
delta_r = reshape(dr, 2, length(dr)/2);

cStruct(2) = Struct;
for ii=1:size(delta_r,2)
    cStruct(2).Vdat(three_verts(ii)).vertxcoord = Struct.Vdat(three_verts(ii)).vertxcoord + delta_r(1,ii);
    cStruct(2).Vdat(three_verts(ii)).vertycoord = Struct.Vdat(three_verts(ii)).vertycoord + delta_r(2,ii);
end

end

