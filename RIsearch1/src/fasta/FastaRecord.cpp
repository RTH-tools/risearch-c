#include "FastaRecord.h"
#include "fasta/fasta.h"


bool FastaRecord::read(FASTAFILE* file)
{
    if (!file || !ReadFASTA(file, m_sequence, m_name)) {
        m_sequence.clear();
        m_name.clear();
        return false;
    }
    return true;
}

const char* FastaRecord::get_sequence() const
{
    return m_sequence.c_str();
}

const char* FastaRecord::get_name() const
{
    return m_name.c_str();
}

std::uint32_t FastaRecord::get_size() const
{
    return static_cast<std::uint32_t>(m_sequence.size());
}
