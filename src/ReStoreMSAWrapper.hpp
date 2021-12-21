#pragma once

// This is a research code proof of concept. I therefore chose to write a header
// only class for simplicity.
#include <cassert>
#include <cstdint>
#include <optional>
#include <vector>

#include <restore/block_serialization.hpp>
#include <restore/core.hpp>
#include <restore/restore_vector.hpp>

#include "MSA.hpp"
#include "PartitionedMSA.hpp"
#include "loadbalance/PartitionAssignment.hpp"

class ReStoreMSAWrapper {
public:
    struct BlockProxy {
        // The information about the current partition is contained in the MSA as each partition has it's own MSA
        // object.
        const MSA* msa;
        uint64_t   siteId;
    };

    ReStoreMSAWrapper(
        const PartitionedMSA& parted_msa, MPI_Comm comm, std::uint16_t replicationLevel,
        const PartitionAssignmentList& partitionToProcessorAssignments) {
        // Parse the partition to processor assignment to get information like the number of ReStore blocks needed etc.
        // Following that, we can create the ReStore object which will hold the MSA data.
        // LOG_DEBUG << "Parsing local processor assignment" << std::endl;
        assert(partitionToProcessorAssignments.size() == ParallelContext::num_ranks());
        for (auto&& assignment: partitionToProcessorAssignments) {
            assert(assignment.length() > 0);
            assert(assignment.begin() != assignment.end());
        }
        // std::cout << "#### Original MSA at rank " << ParallelContext::rank_id() << std::endl;
        // for (auto&& partition: parted_msa.part_list()) {
        //     std::string weights;
        //     for (auto weight: partition.msa().weights()) {
        //         weights += std::to_string(weight) + " ";
        //     }
        //     std::cout << partition.msa().weights().size() << " weights: " << weights << std::endl;
        //     for (auto&& taxon: partition.msa()) {
        //         std::cout << taxon.substr() << std::endl;
        //     }
        // }

        _localPartitionAssignment =
            _parsePartitionToProcessorAssignment(partitionToProcessorAssignments, parted_msa.part_count());
        // Each block consists of one character per taxon plus th weight of this site.
        _restore.emplace(
            comm, replicationLevel, ReStore::OffsetMode::constant, parted_msa.taxon_count() + sizeof(WeightType));

        // Initialize the mapper between ReStore block IDs and <partitions id, site id>
        // LOG_DEBUG << "initialize ReStore block id mapper" << std::endl;
        _restoreBlockMapper.emplace(&_localPartitionAssignment, _partitionIdOffset);
        // LOG_DEBUG << "done initialize ReStore block id mapper" << std::endl;

        // Function to serialize a single site across all taxa as one ReStore block. This sites weight will also be
        // serialized.
        auto serializeBlock = [](const BlockProxy& blockProxy, ReStore::SerializedBlockStoreStream& stream) {
            const auto& msa     = blockProxy.msa;
            const auto& siteId  = blockProxy.siteId;
            const auto& numTaxa = msa->size();

            for (auto taxon = 0u; taxon < numTaxa; ++taxon) {
                // LOG_DEBUG << "Serializing taxon " << taxon << std::endl;
                stream << (*msa)[taxon][siteId];
            }
            stream << (*msa).weights()[siteId];
        };
        // LOG_DEBUG << "Done defining block serializer" << std::endl;

        // Enumerate all Sites in all partitions.
        assert(_localPartitionAssignment.length() > 0);
        assert(_localPartitionAssignment.begin() != _localPartitionAssignment.end());
        auto       currentPartitionRange = _localPartitionAssignment.begin();
        auto       currentSite           = currentPartitionRange->start;
        BlockProxy currentBlockProxy;

        // LOG_DEBUG << "defining site enumerator" << std::endl;
        auto nextBlock = [this, &currentSite, &currentPartitionRange, &parted_msa,
                          endPartitionRange = _localPartitionAssignment.end(),
                          &currentBlockProxy]() -> std::optional<ReStore::NextBlock<BlockProxy>> {
            auto partitionId         = currentPartitionRange->part_id;
            auto currentPartitionEnd = currentPartitionRange->start + currentPartitionRange->length;

            // If we are past the last site of this partition range, advance to the
            // next partition range.
            if (currentSite == currentPartitionEnd) {
                currentPartitionRange++;
                if (currentPartitionRange == endPartitionRange) {
                    return std::nullopt;
                } else {
                    partitionId         = currentPartitionRange->part_id;
                    currentPartitionEnd = currentPartitionRange->length;
                    currentSite         = currentPartitionRange->start;
                }
            }
            assert(currentPartitionRange->length > 0);
            assert(currentSite >= currentPartitionRange->start);
            assert(currentSite < currentPartitionRange->start + currentPartitionRange->length);

            // Compute the ReStore block id for this site.
            const auto restoreBlockId =
                _restoreBlockMapper->partitionAndGlobalSiteId2restoreBlockId(partitionId, currentSite);
            // LOG_DEBUG << "Partition: " << partitionId << " (sites " << currentPartitionRange->start << "-"
            //           << currentPartitionEnd << ") global site: " << currentSite << " ReStore ID: " << restoreBlockId
            //           << std::endl;

            // Assume that the partition IDs between the partition range assignments
            // and the partitions in the Partitioned MSA match up.
            // LOG_DEBUG << "num partitions in MSA: " << parted_msa.part_count()
            //           << "; num partitions in assignment: " << _localPartitionAssignment.num_parts()
            //           << "; num partitions global: " << _numPartitionsGlobal << std::endl;
            assert(_numPartitionsGlobal == parted_msa.part_count());
            assert(currentPartitionRange->part_id < parted_msa.part_count());

            // Get the MSA for this partition.
            const MSA& msa = parted_msa.part_info(partitionId).msa();

            // Assemble return information information, advance to next site and return.
            currentBlockProxy.msa    = &msa;
            currentBlockProxy.siteId = currentSite;
            ReStore::NextBlock<BlockProxy> nextBlock(restoreBlockId, currentBlockProxy);
            currentSite++;
            return nextBlock;
        };

        // LOG_DEBUG << "Start submitting blocks" << std::endl;
        _restore->submitBlocks(serializeBlock, nextBlock, _numReStoreBlocks);
        // LOG_DEBUG << "Done submitting blocks" << std::endl;
    };

    // Delete the copy and move constructors and assignment operators. Copying an MSA is not a good idea, as it uses a
    // lot of memory. Moving might be fine, but I dont need it for now and procastinate thinking about a possible
    // implementation.
    ReStoreMSAWrapper(const ReStoreMSAWrapper& other) = delete;
    ReStoreMSAWrapper& operator=(const ReStoreMSAWrapper& other) = delete;
    ReStoreMSAWrapper(ReStoreMSAWrapper&& other)                 = delete;
    ReStoreMSAWrapper& operator=(ReStoreMSAWrapper&& other) = delete;

    // Update the communicator used by ReStore after a rank failure.
    void updateComm(MPI_Comm newComm) {
        _restore->updateComm(newComm);
    }

    // Populate the given Partitioned MSA with the requested sequences. Use the
    // data from the ReStore. Uses the local sequence ranges stored in the MSA
    // object. (For example the RangeList passed to the constructor.) It is only
    // valid to pass MSA object which do not contain any sequences.
    void restoreMSA(PartitionedMSA& parted_msa, const PartitionAssignmentList& partitionToProcessorAssignments) {
        // Update the local partition assignment. As the _restoreBlockMapper has a reference to our
        // _localPartitionAssignment, we do not need to explicitly update the mapper.
        _localPartitionAssignment = partitionToProcessorAssignments[ParallelContext::rank_id()];
        assert(_localPartitionAssignment.length() > 0);
        _restoreBlockMapper->updateThisRanksPartitionAssignment(&_localPartitionAssignment);

        // Build the block request data structure. It defines which ReStore blocks
        // need to be moved to which rank.
        using BlockRange   = std::pair<ReStore::block_id_t, size_t>;
        using BlockRequest = std::pair<BlockRange, ReStoreMPI::current_rank_t>;
        std::vector<BlockRequest> blockRequests;

        auto numRanks = asserting_cast<ReStoreMPI::current_rank_t>(ParallelContext::num_ranks());
        for (ReStoreMPI::current_rank_t destRank = 0; destRank < numRanks; destRank++) {
            const auto& partitionAssignment = partitionToProcessorAssignments.at(destRank);
            for (auto& partition: partitionAssignment) {
                auto restoreBlockIdOfFirstSite =
                    _restoreBlockMapper->partitionAndGlobalSiteId2restoreBlockId(partition.part_id, partition.start);

                blockRequests.emplace_back(BlockRange(restoreBlockIdOfFirstSite, partition.length), destRank);
                // LOG_DEBUG << "Requesting blocks " << restoreBlockIdOfFirstSite << "-"
                //           << restoreBlockIdOfFirstSite + partition.length - 1 << " to rank "
                //           << destRank << std::endl;
            }
        }

        // Define the callback function ReStore uses to hand us the received block in serialized form. This function
        // deserialized the incoming blocks and places the contained sites into a local store for the sequences. In a
        // second step, we then populate the MSA object. This has to be done, as writing single sites to a partitioned
        // MSA is not trivial and not currently supported by the PartitionedMSA object.
        // TODO How to get the number of taxa to reserve enough space beforehand?
        // TODO Does this work?
        auto numTaxa = parted_msa.part_list().at(0).msa().size();
        assert(numTaxa > 0);

        // Caluclate the number of local sites per partition and reserve space for the received sequences.
        std::vector<size_t>                         numLocalSitesPerPartition(_numPartitionsGlobal, 0);
        std::vector<std::vector<std::vector<char>>> sequences(_numPartitionsGlobal); // By partition ID and taxon
        std::vector<WeightVector>                   weights(_numPartitionsGlobal);   // But partition ID
        assert(numLocalSitesPerPartition.size() == _numPartitionsGlobal);
        assert(sequences.size() == _numPartitionsGlobal);
        assert(weights.size() == sequences.size());

        assert(_localPartitionAssignment.length() > 0);
        for (auto& partition: _localPartitionAssignment) {
            // Number of sites per partition
            auto partitionId     = partition.part_id;
            auto partitionLength = partition.length;
            assert(partitionLength > 0);
            assert(numLocalSitesPerPartition.at(partitionId) == 0);
            // LOG_DEBUG << "Partition " << partitionId << " has " << partitionLength << " local sites" << std::endl;
            numLocalSitesPerPartition[partitionId] = partitionLength;

            // Reserve space for the received sequences. Actually resize() instead of reserve(), as we are doing
            // assignment during the block deserialization.
            sequences[partitionId].resize(numTaxa);
            weights[partitionId].resize(numLocalSitesPerPartition[partitionId]);
            for (auto taxon = 0u; taxon < numTaxa; ++taxon) {
                sequences[partitionId][taxon].resize(numLocalSitesPerPartition[partitionId]);
            }
        }

        auto blockDeserializer = [this, numTaxa, &weights,
                                  &sequences](const void* data, size_t lengthInBytes, ReStore::block_id_t blockId) {
            // LOG_DEBUG << "Received block " << blockId << " at rank " << ParallelContext::rank_id() << std::endl;
            assert(_restoreBlockMapper);
            auto [partitionId, localSiteId] = _restoreBlockMapper->restoreBlockId2PartitionAndLocalSiteId(blockId);
            // LOG_DEBUG << "Block belongs to partition " << partitionId << " and local site " << localSiteId <<
            // std::endl;

            // Copy the characters to their respective taxon.
            assert(lengthInBytes == numTaxa * sizeof(char) + sizeof(WeightType));
            for (auto taxon = 0u; taxon < numTaxa; ++taxon) {
                // LOG_DEBUG << "Copying site " << localSiteId << " of taxon " << taxon << " to partition " <<
                // partitionId
                //           << ": " << *(reinterpret_cast<const char*>(data) + taxon) << std::endl;
                assert(sequences.size() > partitionId);
                assert(taxon < sequences[partitionId].size());
                assert(sequences[partitionId][taxon].size() > 0);
                assert(localSiteId < sequences[partitionId][taxon].size());
                sequences[partitionId][taxon][localSiteId] = *(reinterpret_cast<const char*>(data) + taxon);
            }

            // Copy over the sites' weight.
            // I want to do the pointer arithmetic on a char pointer.
            const WeightType* weightPtr =
                reinterpret_cast<const WeightType*>(reinterpret_cast<const char*>(data) + numTaxa);
            weights[partitionId][localSiteId] = *weightPtr;
        };

        // for (auto& blockRequest: blockRequests) {
        //     auto& [blockRange, destRank] = blockRequest;
        //     LOG_DEBUG << "Requesting blocks " << blockRange.first << "-" << blockRange.first + blockRange.second - 1
        //               << " to rank " << destRank << std::endl;
        // }
        // Fetch the MSA sites from the ReStore into the sequences store.
        _restore->pushBlocksCurrentRankIds(blockRequests, blockDeserializer);

        // Populate the MSA object with the sequences.
        for (size_t partitionId = 0; partitionId < _numPartitionsGlobal; partitionId++) {
            auto      partitionAssignment = _localPartitionAssignment.find(partitionId);
            RangeList rangeList;
            rangeList.emplace_back(partitionAssignment->start, partitionAssignment->length);
            auto& msa = parted_msa.part_list().at(partitionId).msa() = MSA(rangeList);
            // LOG_DEBUG << "Populating MSA for partition " << partitionId << " with "
            //           << numLocalSitesPerPartition[partitionId] << " sites and " << sequences[partitionId].size()
            //           << " taxa." << std::endl;
            for (auto& sequenceOfTaxon: sequences[partitionId]) {
                sequenceOfTaxon.push_back('\0');
                // LOG_DEBUG << "Adding sequence of length " << sequenceOfTaxon.size() - 1 << " to partition "
                //           << partitionId << std::endl;
                // LOG_DEBUG << std::string(sequenceOfTaxon.data()) << std::endl;
                msa.append(std::string(sequenceOfTaxon.data()));
            }
            msa.weights(std::move(weights[partitionId]));
        }
        // std::cout << "#### Restored MSA at rank " << ParallelContext::rank_id() << std::endl;
        // for (auto&& partition: parted_msa.part_list()) {
        //     std::string weights;
        //     for (auto weight: partition.msa().weights()) {
        //         weights += std::to_string(weight) + " ";
        //     }
        //     std::cout << partition.msa().weights().size() << " weights: " << weights << std::endl;
        //     for (auto&& taxon: partition.msa()) {
        //         std::cout << taxon.substr() << std::endl;
        //     }
        // }
    }

private:
    // This class encapsulates the functionality used to map partition ids + site ids to ReStore ids and vice versa. It
    // also maps global to local site ids and vice versa. Global site ids are a site's position in the partition as in
    // the MSA. Local site ids are the position of this site in this ranks MSA object.
    class ReStoreBlockMapper {
        // The global site id is this sites unique index across all sites in the whole MSA instead of only on this rank.
    public:
        ReStoreBlockMapper(
            const PartitionAssignment* thisRanksPartitionAssignment, const std::vector<size_t>& partitionIdOffsets)
            : _partitionAssignment(thisRanksPartitionAssignment),
              _partitionIdOffsets(partitionIdOffsets) {}

        void updateThisRanksPartitionAssignment(const PartitionAssignment* thisRanksPartitionAssignment) {
            _partitionAssignment = thisRanksPartitionAssignment;
        }

        inline std::pair<uint64_t, uint64_t> restoreBlockId2PartitionAndGlobalSiteId(uint64_t restoreBlockId) {
            // Determine to which partition this block id belongs.
            auto partitionId =
                std::distance(
                    _partitionIdOffsets.begin(),
                    std::upper_bound(_partitionIdOffsets.begin(), _partitionIdOffsets.end(), restoreBlockId))
                - 1;

            auto globalSiteIdInPartition = restoreBlockId - _partitionIdOffsets[partitionId];

            return std::pair(partitionId, globalSiteIdInPartition);
        }

        inline uint64_t partitionAndGlobalSiteId2restoreBlockId(uint64_t partitionId, uint64_t globalSiteId) {
            assert(partitionId < _partitionIdOffsets.size());
            return _partitionIdOffsets[partitionId] + globalSiteId;
        }

        inline uint64_t partitionAndLocalSiteId2restoreBlockId(uint64_t partitionId, uint64_t localSiteId) {
            assert(partitionId < _partitionIdOffsets.size());
            return _partitionIdOffsets[partitionId] + local2globalSiteId(partitionId, localSiteId);
        }

        inline std::pair<uint64_t, uint64_t> restoreBlockId2PartitionAndLocalSiteId(uint64_t restoreBlockId) {
            auto result   = restoreBlockId2PartitionAndGlobalSiteId(restoreBlockId);
            result.second = global2localSiteId(result.first, result.second);
            return result;
        }

        inline uint64_t local2globalSiteId(uint64_t partitionId, uint64_t localSiteId) {
            assert(_partitionAssignment->length() > 0);
            auto partition = _partitionAssignment->find(partitionId);
            assert(partition != _partitionAssignment->end());
            return partition->start + localSiteId;
        }

        inline uint64_t global2localSiteId(uint64_t partitionId, uint64_t globalSiteId) {
            assert(_partitionAssignment->length() > 0);
            auto partition = _partitionAssignment->find(partitionId);
            // for (auto& partition: *_partitionAssignment) {
            //     LOG_DEBUG << "My partition assignment: Partition " << partition.part_id << ": start at "
            //               << partition.start << " with length " << partition.length << std::endl;
            // }
            assert(partition != _partitionAssignment->end());
            return globalSiteId - partition->start;
        }

    private:
        const PartitionAssignment* _partitionAssignment; // Needs to be updated after rank failure or rebalancing.
        const std::vector<size_t>& _partitionIdOffsets;
    };

    std::optional<ReStoreBlockMapper> _restoreBlockMapper;

    // The ReStore object we use to store the sequencing data.
    std::optional<ReStore::ReStore<BlockProxy>> _restore;

    // The number of partitions in our analysis.
    uint64_t _numPartitionsGlobal;
    uint64_t _numPartitionsLocal; // on this rank

    // The number of ReStore blocks. In the current implementation, this is the number of sites across all partitions.
    uint64_t _numReStoreBlocks;

    // A sites (ReStore) block id is computed by its partition id offset plus it's position in the given partition.
    // _partIdOffset therefore stores the exclusive prefix sum of the partition sizes.
    std::vector<uint64_t> _partitionIdOffset;

    // The partition ranges assigned to this rank.
    PartitionAssignment _localPartitionAssignment;

    // Parses the partition assignment list to extract the following information:
    // - The number of partitions
    // - The number of ReStore blocks
    // - The partitions id offsets
    // - The local partition assignment -> return value
    const PartitionAssignment& _parsePartitionToProcessorAssignment(
        const PartitionAssignmentList& partitionToProcessorAssignments, size_t numPartitionsInMSA) {
        // We are working in MPI only mode (one thread per rank).;
        assert(partitionToProcessorAssignments.size() == ParallelContext::num_ranks());
        _numPartitionsGlobal = 0;
        _numPartitionsLocal  = 0;
        _numReStoreBlocks    = 0;

        size_t numSites = 0;
        _partitionIdOffset.resize(numPartitionsInMSA, 0);
        std::vector<size_t> partitionSizes(numPartitionsInMSA, 0);
        // LOG_DEBUG << "There are " << numPartitionsInMSA << " partitions in the MSA." << std::endl;

        for (auto partAssignment: partitionToProcessorAssignments) {
            for (auto partRange: partAssignment) {
                // First, collect the sizes of all partitions in partitionSizes
                // Only later compute the exclusive prefix sum.
                assert(partRange.part_id < partitionSizes.size());
                partitionSizes[partRange.part_id] += partRange.length;

                // Assume, that the partition ids are consecutive.
                _numPartitionsGlobal = std::max(_numPartitionsGlobal, partRange.part_id);

                // Also count the number of sites to later infer the number of ReStore
                // blocks. Assume, that there are no overlapping partitions.
                numSites += partRange.length;
            }
        }
        _numPartitionsGlobal++;
        // LOG_DEBUG << "There are " << _numPartitionsGlobal << " global partitions." << std::endl;
        assert(_numPartitionsGlobal == _partitionIdOffset.size());
        assert(_numPartitionsGlobal == partitionSizes.size());

        // Compute the exclusive prefix sum over the partition sizes to obtain the
        // partition's id offset.
        assert(_partitionIdOffset.size() == partitionSizes.size());
        std::exclusive_scan(partitionSizes.begin(), partitionSizes.end(), _partitionIdOffset.begin(), 0);
        assert(_partitionIdOffset.size() == _numPartitionsGlobal);

        // for (size_t partitionId = 0; partitionId < _numPartitionsGlobal; partitionId++) {
        //     LOG_DEBUG << "Partition " << partitionId << ": start at " << _partitionIdOffset[partitionId]
        //               << " with length " << partitionSizes[partitionId] << std::endl;
        // }

        // Currently, the number of ReStore blocks is equal to the number of sites across all partitions in the MSA.
        _numReStoreBlocks = numSites;

        const auto& localPartitionAssignment = partitionToProcessorAssignments[ParallelContext::rank_id()];
        _numPartitionsLocal                  = localPartitionAssignment.num_parts();
        assert(_numPartitionsLocal > 0);
        assert(_numPartitionsLocal <= _numPartitionsGlobal);
        return localPartitionAssignment;
    }
};