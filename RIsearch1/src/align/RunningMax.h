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

    void set_if_better(int candidate, int pos_i, int pos_j)
    {
        if (candidate > score) {
            score = candidate;
            this->pos_i = pos_i;
            this->pos_j = pos_j;
        }
    }
};

struct RunningRowMax {
    int score;
    int pos_i;

    void set(int candidate, int pos_col)
    {
        score = candidate;
        pos_i = pos_col;
    }

    void set_if_better(int candidate, int pos_col)
    {
        if (candidate > score) {
            score = candidate;
            pos_i = pos_col;
        }
    }
};
