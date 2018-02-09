N = 10:24;
nIter = 20;
for ii = 1:length(N)
    for n = 1:nIter
%         N = 13;
        alpha = .02;
        beta = 2;

        [ q, ~, b0 ] = generate.meshGrid(N(ii),alpha);
        q(b0,:) = 1.05*q(b0,:);
        p = 3*ones(size(q,1),1);
        z = beta*( rand(size(p)) - .5 ); 
        p = p + z;
        % p(b0) = 1;
        theta = .1 * ( rand(size(p))) ;

        [ tri, L, ~, q, p, theta, contin ] = pressure.calculateTri( q, p, theta, b0 );
        [L,Struct] = seg.generate_structs(L,0,1);
        Struct = seg.threefold_cell(Struct);
        Struct = seg.recordBonds(Struct,L);
        Struct = seg.curvature(Struct,size(L));

        Struct.labelMat = L;
        tic
        [PN,ERes,r0] = fitDual.returnDual(Struct,3,1);
        delta(ii,n) = toc;
        numCells(ii,n) = size(PN{1}.q,1);
    end
end
% clearvars -except delta numCells
% PN = pressure.net(q,theta,p,tri);
% kStruct = PN.returnStruct();
% kStruct = seg.threefold_cell(kStruct);

% clearvars -except L kStruct Struct PN
% clear z b0
% 
% % etaVec = linspace(0,.1,10);
% etaVec = 0;
% n = 1;
% for eta = etaVec
%     Struct = PN.returnStruct(eta);
%     Struct = seg.threefold_cell(Struct);
% 
%     oStruct = MI.storeMech(Struct,1,3);
%     % PNi = fitDual.returnDual(Struct,3,1);
%     % nStruct = PNi{1}.uploadMechanics(Struct);
% 
%     To = [];
%     aTo = [];
%     for b = 1:length(oStruct.Bdat)
%         if (~isempty(oStruct.Bdat(b).tension) && ~isempty(oStruct.Bdat(b).actual_tension) )
%             To = [To,oStruct.Bdat(b).tension];
%             aTo = [aTo,oStruct.Bdat(b).actual_tension];
%         end
%     end
% 
%     Po = [];
%     aPo = [];
%     for c = 1:length(oStruct.Cdat)
%         if (~isempty(oStruct.Cdat(c).pressure) && ~isempty(oStruct.Cdat(c).actual_pressure) )
%             Po = [Po,oStruct.Cdat(c).pressure];
%             aPo = [aPo,oStruct.Cdat(c).actual_pressure];
%         end
%     end
%     
%     cT(n) = corr(To',aTo');
%     cP(n) = corr(Po',aPo');
%     n = n + 1;
% end
% 
% % Tn = [];
% % aTn = [];
% % for b = 1:length(nStruct.Bdat)
% %     if (~isempty(nStruct.Bdat(b).tension) && ~isempty(nStruct.Bdat(b).actual_tension) )
% %         Tn = [Tn,nStruct.Bdat(b).tension];
% %         aTn = [aTn,nStruct.Bdat(b).actual_tension];
% %     end
% % end