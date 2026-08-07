#pragma once

static void set_alignment_symbols (char query_nt, char target_nt, char *query_alignment,
               char *target_alignment)
{
    /*this function takes two nucleotides as inputs and set characters in the alignment string to:
     | : if the two nt can base pair,
     M : if there is a mismatch,
     W : if there is a wobble base pair*/
    if ('A' == query_nt || 'A' == target_nt)
    {
        if ('U' == query_nt || 'U' == target_nt)
        {
            *query_alignment = '|';
            *target_alignment = '|';
        }
        else
        {
            *query_alignment = 'M';
            *target_alignment = 'M';
        }
    }
    else if ('G' == query_nt || 'G' == target_nt)
    {
        if ('U' == query_nt || 'U' == target_nt)
        {
            *query_alignment = 'W';
            *target_alignment = 'W';
        }
        else if ('C' == query_nt || 'C' == target_nt)
        {
            *query_alignment = '|';
            *target_alignment = '|';
        }
        else
        {
            *query_alignment = 'M';
            *target_alignment = 'M';
        }
    }
    else
    {
        *query_alignment = 'M';
        *target_alignment = 'M';
    }

}
