#include "FastaRecord.h"

const unsigned char* FastaRecord::get_sequence() const
{
    return m_sequence.get();
}

const unsigned char* FastaRecord::get_name() const
{
    return m_name.get();
}
unsigned int FastaRecord::get_size() const
{
    return m_size;
}