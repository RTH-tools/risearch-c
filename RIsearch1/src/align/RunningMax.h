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
