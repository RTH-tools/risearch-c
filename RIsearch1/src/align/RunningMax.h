#pragma once


struct RunningMax {
    int score;
    int pos_i;
    int pos_j;

    void set(int candidate, int pos_i, int pos_j)
    {
        score = candidate;
        this->pos_i = pos_i;
        this->pos_j = pos_j;
    }

    /* Whether candidate would displace what is held. Strictly better, so the
       first of equal candidates keeps the position. */
    bool test(int candidate) const
    {
        return candidate > score;
    }

    void set_if_better(int candidate, int pos_i, int pos_j)
    {
        if (test(candidate)) {
            set(candidate, pos_i, pos_j);
        }
    }
};

struct RunningVectorMax {
    int score;
    int pos_i;

    void set(int candidate, int pos_col)
    {
        score = candidate;
        pos_i = pos_col;
    }

    bool test(int candidate) const
    {
        return candidate > score;
    }

    void set_if_better(int candidate, int pos_col)
    {
        if (test(candidate)) {
            set(candidate, pos_col);
        }
    }
};
