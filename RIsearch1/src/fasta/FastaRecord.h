#pragma once

#include "memory/MallocRAII.hpp"

class FastaRecord {
public:
    [[nodiscard]] const unsigned char* get_sequence() const;
    [[nodiscard]] const unsigned char* get_name() const;
    [[nodiscard]] unsigned int get_size() const;


private:
    MallocRAII<unsigned char> m_sequence;
    MallocRAII<unsigned char> m_name;
    unsigned int m_size = 0;
};
