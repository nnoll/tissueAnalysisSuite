L = seg.memWS(mem,50,0,1,3.5);
[L,Struct] = seg.generate_structs(L,0,1);
L = seg.removeBadCells(Struct,L);
L = seg.relabelL(L);
[L,Struct] = seg.generate_structs(L,0,0);

Struct = seg.threefold_cell(Struct);
Struct = seg.recordBonds(Struct,L);
Struct = seg.curvature(Struct,size(L));
Struct = seg.removeFourFold(Struct,size(L));

Struct = seg.makeConvexArray(Struct);
cStruct = isogonal.imposeComptCond(Struct,.5);
cStruct = seg.makeConvexArray(cStruct);


