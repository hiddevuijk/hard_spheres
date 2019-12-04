#ifndef GUARD_CELL_INDEX_H
#define GUARD_CELL_INDEX_H


class Cindex {
public:
    Cindex(): i(-1), j(-1), k(-1) {}
    Cindex(int ii, int jj, int kk): i(ii), j(jj), k(kk) {}
    int i,j,k;
};

#endif
