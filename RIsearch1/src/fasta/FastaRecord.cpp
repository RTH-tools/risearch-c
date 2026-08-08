#include "FastaRecord.h"
#include "fasta/fasta.h"


bool FastaRecord::read(FASTAFILE* file)
{
    char* name = nullptr;
    char* sequence = nullptr;
    std::uint32_t size = 0;

    if (!file || !ReadFASTA(file, &sequence, &name, &size)) {
        m_name.reset();
        m_sequence.reset();
        m_size = 0;
        return false;
    }
    m_name.reset(name);
    m_sequence.reset(sequence);
    m_size = size;
    return true;
}

const char* FastaRecord::get_sequence() const
{
    return m_sequence.get();
}

const char* FastaRecord::get_name() const
{
    return m_name.get();
}

std::uint32_t FastaRecord::get_size() const
{
    return m_size;
}