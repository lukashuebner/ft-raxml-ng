

#include "RaxmlTest.hpp"
#include "src/ReStoreMSAWrapper.hpp"
//#include "PartitioneDMSA.hpp"

using namespace std;

/*
TEST(ReStoreMSAIntegration, OnePartition) {
    constexpr uint64_t numPartitions = 1;

    // Populate the MSA object with the sequences.
    PartitionedMSA parted_msa;

    PartitionInfo partition1_info;
    RangeList partition1_range_list;
    partition1_range_list.emplace_back(0, 4);
    MSA partition1_msa(partition1_range_list);
    partition1_msa.append("AAAA");
    partition1_msa.append("CCCC");
    partition1_info.msa(std::move(partition1_msa));
    
    // Build the partition assignment list

    // Create the ReStore object
    ReStoreMSAWrapper restoreWrapper(parted_msa, MPI_COMM_WORLD, 2, PartitionAssignment);
}

TEST(ReStoreMSAIntegration, OnePartition) {
    constexpr uint64_t numPartitions = 1;

    // Populate the MSA object with the sequences.
    for (size_t partitionId = 0; partitionId < _numPartitionsGlobal; partitionId++) {
        auto      partitionAssignment = _localPartitionAssignment.find(partitionId);
        RangeList rangeList;
        rangeList.emplace_back(partitionAssignment->start, partitionAssignment->length);
        auto& msa = parted_msa.part_list().at(partitionId).msa() = MSA(rangeList);
        LOG_DEBUG << "Populating MSA for partition " << partitionId << " with "
                  << numLocalSitesPerPartition[partitionId] << " sites and " << sequences[partitionId].size()
                  << " taxa." << std::endl;
        for (auto& sequenceOfTaxon: sequences[partitionId]) {
            sequenceOfTaxon.push_back('\0');
            LOG_DEBUG << "Adding sequence of length " << sequenceOfTaxon.size() - 1 << " to partition " << partitionId
                      << std::endl;
            LOG_DEBUG << std::string(sequenceOfTaxon.data()) << std::endl;
            msa.append(std::string(sequenceOfTaxon.data()));
        }
    }
    std::cout << "#### Restored MSA at rank " << ParallelContext::rank_id() << std::endl;
    for (auto&& partition: parted_msa.part_list()) {
        for (auto&& taxon: partition.msa()) {
            std::cout << taxon.substr() << std::endl;
        }
    }
}
*/