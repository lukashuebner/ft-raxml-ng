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
        const MSA* msaOfThisPattern;
        uint64_t   localPatternId;
    };

    ReStoreMSAWrapper(
        const PartitionedMSA& parted_msa, MPI_Comm comm, std::uint16_t replicationLevel,
        const PartitionAssignmentList& partitionToProcessorAssignments) {
        // Parse the partition to processor assignment to get information like the number of ReStore blocks needed etc.
        // Following that, we can create the ReStore object which will hold the MSA data.
        assert(partitionToProcessorAssignments.size() == ParallelContext::num_ranks());
        for (auto&& assignment: partitionToProcessorAssignments) {
            assert(assignment.length() > 0);
            assert(assignment.begin() != assignment.end());
        }

        _localPartitionAssignment =
            _parsePartitionToProcessorAssignment(partitionToProcessorAssignments, parted_msa.part_count());
        // Each block consists of one character per taxon plus th weight of this pattern.
        _restore.emplace(
            comm, replicationLevel, ReStore::OffsetMode::constant, parted_msa.taxon_count() + sizeof(WeightType));

        // Initialize the mapper between ReStore block IDs and <partitions id, pattern id>
        _restoreBlockMapper.emplace(&_localPartitionAssignment, _partitionIdOffset);

        // Function to serialize a single pattern across all taxa as one ReStore block. This pattern's weight will also
        // be serialized.
        auto serializeBlock = [](const BlockProxy& blockProxy, ReStore::SerializedBlockStoreStream& stream) {
            const auto& msa       = blockProxy.msaOfThisPattern;
            const auto& patternId = blockProxy.localPatternId;
            const auto& numTaxa   = msa->size();

            for (auto taxon = 0u; taxon < numTaxa; ++taxon) {
                assert((*msa)[taxon][patternId] != '\0');
                stream << (*msa)[taxon][patternId];
            }
            assert(patternId < msa->weights().size());
            assert(msa->weights()[patternId] != 0);
            stream << msa->weights()[patternId];
        };

        // Enumerate all patterns in all partitions.
        assert(_localPartitionAssignment.length() > 0);
        assert(_localPartitionAssignment.begin() != _localPartitionAssignment.end());
        auto       partitionRange = _localPartitionAssignment.begin();
        auto       partitionId    = partitionRange->part_id;
        auto       partitionEnd   = partitionRange->start + partitionRange->length;
        auto       patternId      = partitionRange->start;
        const MSA* msa            = &(parted_msa.part_info(partitionId).msa());
        BlockProxy blockProxy;

        auto nextBlock = [this, &patternId, &partitionRange, &parted_msa, &msa, &partitionEnd, &partitionId,
                          localPartitionAssignmentEnd = _localPartitionAssignment.end(),
                          &blockProxy]() -> std::optional<ReStore::NextBlock<BlockProxy>> {
            // If we are past the last pattern of this partition range, advance to the next partition range.
            if (patternId == partitionEnd) {
                partitionRange++;
                if (partitionRange == localPartitionAssignmentEnd) {
                    return std::nullopt;
                } else {
                    partitionId  = partitionRange->part_id;
                    partitionEnd = partitionRange->start + partitionRange->length;
                    patternId    = partitionRange->start;
                    msa          = &(parted_msa.part_info(partitionId).msa());
                }
            }
            assert(partitionRange->length > 0);
            assert(patternId >= partitionRange->start);
            assert(patternId < partitionEnd);

            // Compute the ReStore block id for this pattern.
            const auto restoreBlockId =
                _restoreBlockMapper->partitionAndGlobalPatternId2restoreBlockId(partitionId, patternId);

            // Assume that the partition IDs between the partition range assignments
            // and the partitions in the Partitioned MSA match up.
            assert(_numPartitionsGlobal == parted_msa.part_count());
            assert(partitionRange->part_id < parted_msa.part_count());

            // Assemble return information information, advance to next pattern and return.
            blockProxy.msaOfThisPattern = msa;
            blockProxy.localPatternId   = _restoreBlockMapper->global2localPatternId(partitionId, patternId);
            ReStore::NextBlock<BlockProxy> nextBlock(restoreBlockId, blockProxy);
            patternId++;
            return nextBlock;
        };

        _restore->submitBlocks(serializeBlock, nextBlock, _numReStoreBlocks);
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

        auto numRanks           = asserting_cast<ReStoreMPI::current_rank_t>(ParallelContext::num_ranks());
        auto numBlocksRequested = 0u;
        for (ReStoreMPI::current_rank_t destRank = 0; destRank < numRanks; destRank++) {
            const auto& partitionAssignment = partitionToProcessorAssignments.at(destRank);
            for (auto& partition: partitionAssignment) {
                auto restoreBlockIdOfFirstPattern =
                    _restoreBlockMapper->partitionAndGlobalPatternId2restoreBlockId(partition.part_id, partition.start);

                blockRequests.emplace_back(BlockRange(restoreBlockIdOfFirstPattern, partition.length), destRank);

                if (destRank == static_cast<int>(ParallelContext::rank_id())) {
                    numBlocksRequested += partition.length;
                }
            }
        }

        // Define the callback function ReStore uses to hand us the received block in serialized form. This function
        // deserialized the incoming blocks and places the contained patterns into a local store for the sequences. In a
        // second step, we then populate the MSA object. This has to be done, as writing single patterns to a
        // partitioned MSA is not trivial and not currently supported by the PartitionedMSA object.
        auto numTaxa = parted_msa.taxon_count();
        assert(numTaxa > 0);

        // Caluclate the number of local patterns per partition and reserve space for the received sequences.
        std::vector<size_t>                         numLocalPatternsPerPartition(_numPartitionsGlobal, 0);
        std::vector<std::vector<std::vector<char>>> sequences(_numPartitionsGlobal); // By partition ID and taxon
        std::vector<WeightVector>                   weights(_numPartitionsGlobal);   // But partition ID
        assert(numLocalPatternsPerPartition.size() == _numPartitionsGlobal);
        assert(sequences.size() == _numPartitionsGlobal);
        assert(weights.size() == sequences.size());

        assert(_localPartitionAssignment.length() > 0);
        auto sizeOfSequenceBuffer = 0u;
        for (auto& partition: _localPartitionAssignment) {
            // Number of patterns per partition
            auto partitionId     = partition.part_id;
            auto partitionLength = partition.length;
            assert(partitionLength > 0);
            assert(numLocalPatternsPerPartition[partitionId] == 0);
            numLocalPatternsPerPartition[partitionId] = partitionLength;
            assert(numLocalPatternsPerPartition[partitionId] > 0);

            // Reserve space for the received sequences. Actually resize() instead of reserve(), as we are doing
            // assignment during the block deserialization.
            sequences[partitionId].resize(numTaxa);
            weights[partitionId].resize(numLocalPatternsPerPartition[partitionId]);
            sizeOfSequenceBuffer += numLocalPatternsPerPartition[partitionId];
            for (auto taxon = 0u; taxon < numTaxa; ++taxon) {
                sequences[partitionId][taxon].resize(numLocalPatternsPerPartition[partitionId]);
            }
        }

        auto numBlocksReceived = 0u;
        auto blockDeserializer = [this, numTaxa, &weights, &numBlocksReceived,
                                  &sequences](const void* data, size_t lengthInBytes, ReStore::block_id_t blockId) {
            numBlocksReceived++;
            assert(_restoreBlockMapper);
            auto [partitionId, localPatternId] =
                _restoreBlockMapper->restoreBlockId2PartitionAndLocalPatternId(blockId);

            // Copy the characters to their respective taxon.
            assert(lengthInBytes == numTaxa * sizeof(char) + sizeof(WeightType));
            for (auto taxon = 0u; taxon < numTaxa; ++taxon) {
                assert(sequences.size() > partitionId);
                assert(taxon < sequences[partitionId].size());
                assert(sequences[partitionId][taxon].size() > 0);
                assert(localPatternId < sequences[partitionId][taxon].size());
                assert(sequences[partitionId][taxon][localPatternId] == '\0'); // No pattern should be written twice.
                sequences[partitionId][taxon][localPatternId] = *(reinterpret_cast<const char*>(data) + taxon);
                assert(sequences[partitionId][taxon][localPatternId] != '\0'); // That's not a valid IUPAC code.
            }

            // Copy over the pattern's weight.
            // I want to do the pointer arithmetic on a char pointer.
            const WeightType* weightPtr =
                reinterpret_cast<const WeightType*>(reinterpret_cast<const char*>(data) + numTaxa);
            assert(*weightPtr > 0);
            assert(localPatternId < weights[partitionId].size());
            weights[partitionId][localPatternId] = *weightPtr;
        };

        // Fetch the MSA patterns from the ReStore into the sequences store.
        _restore->pushBlocksCurrentRankIds(blockRequests, blockDeserializer);

        assert(numBlocksReceived == numBlocksRequested);
        assert(numBlocksReceived == sizeOfSequenceBuffer);

        // Populate the MSA object with the sequences.
        for (size_t partitionId = 0; partitionId < _numPartitionsGlobal; partitionId++) {
            // Get the local partition assignment for this partition. If this rank did not get any patterns of this
            // partition, we do not need to populate our MSA with sequences either and can therefore skip it.
            auto partitionAssignment = _localPartitionAssignment.find(partitionId);
            if (partitionAssignment == _localPartitionAssignment.end()) {
                continue;
            }

            // Build the RangeList object describing our assignment of this partition.
            RangeList rangeList;
            rangeList.emplace_back(partitionAssignment->start, partitionAssignment->length);
            auto& msa = parted_msa.part_list().at(partitionId).msa() = MSA(rangeList);

            // Append the sequences (one for each taxon) of our assignment of this partition to the corresponding MSA
            // object.
            for (auto& sequenceOfTaxon: sequences[partitionId]) {
                sequenceOfTaxon.push_back('\0');
                std::string sequenceOfTaxon_str(sequenceOfTaxon.begin(), sequenceOfTaxon.end());
                assert(sequenceOfTaxon_str.length() == sequenceOfTaxon.size());
                msa.append(sequenceOfTaxon_str);
            }
            assert(partitionId < weights.size());
            msa.weights(std::move(weights[partitionId]));
        }
    }

private:
    // This class encapsulates the functionality used to map partition ids + pattern ids to ReStore ids and vice versa.
    // It also maps global to local pattern ids and vice versa. Global pattern ids are a pattern's position in the
    // partition as in the MSA. Local pattern ids are the position of this pattern in this ranks MSA object.
    class ReStoreBlockMapper {
        // The global pattern id is this pattern's unique index across all patterns in the whole MSA instead of only on
        // this rank.
    public:
        ReStoreBlockMapper(
            const PartitionAssignment* thisRanksPartitionAssignment, const std::vector<size_t>& partitionIdOffsets)
            : _partitionAssignment(thisRanksPartitionAssignment),
              _partitionIdOffsets(partitionIdOffsets) {}

        void updateThisRanksPartitionAssignment(const PartitionAssignment* thisRanksPartitionAssignment) {
            _partitionAssignment = thisRanksPartitionAssignment;
        }

        inline std::pair<uint64_t, uint64_t> restoreBlockId2PartitionAndGlobalPatternId(uint64_t restoreBlockId) {
            // Determine to which partition this block id belongs.
            auto partitionId =
                std::distance(
                    _partitionIdOffsets.begin(),
                    std::upper_bound(_partitionIdOffsets.begin(), _partitionIdOffsets.end(), restoreBlockId))
                - 1;

            auto globalPatternIdInPartition = restoreBlockId - _partitionIdOffsets[partitionId];

            return std::pair(partitionId, globalPatternIdInPartition);
        }

        inline uint64_t partitionAndGlobalPatternId2restoreBlockId(uint64_t partitionId, uint64_t globalPatternId) {
            assert(partitionId < _partitionIdOffsets.size());
            return _partitionIdOffsets[partitionId] + globalPatternId;
        }

        inline uint64_t partitionAndLocalPatternId2restoreBlockId(uint64_t partitionId, uint64_t localPatternId) {
            assert(partitionId < _partitionIdOffsets.size());
            return _partitionIdOffsets[partitionId] + local2globalPatternId(partitionId, localPatternId);
        }

        inline std::pair<uint64_t, uint64_t> restoreBlockId2PartitionAndLocalPatternId(uint64_t restoreBlockId) {
            auto result   = restoreBlockId2PartitionAndGlobalPatternId(restoreBlockId);
            result.second = global2localPatternId(result.first, result.second);
            return result;
        }

        inline uint64_t local2globalPatternId(uint64_t partitionId, uint64_t localPatternId) {
            assert(_partitionAssignment->length() > 0);
            auto partition = _partitionAssignment->find(partitionId);
            assert(partition != _partitionAssignment->end());
            return partition->start + localPatternId;
        }

        inline uint64_t global2localPatternId(uint64_t partitionId, uint64_t globalPatternId) {
            assert(_partitionAssignment->length() > 0);
            auto partition = _partitionAssignment->find(partitionId);
            assert(partition != _partitionAssignment->end());
            return globalPatternId - partition->start;
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

    // The number of ReStore blocks. In the current implementation, this is the number of patterns across all
    // partitions.
    uint64_t _numReStoreBlocks;

    // A pattern's (ReStore) block id is computed by its partition id offset plus it's position in the given partition.
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

        size_t numPatterns = 0;
        _partitionIdOffset.resize(numPartitionsInMSA, 0);
        std::vector<size_t> partitionSizes(numPartitionsInMSA, 0);

        for (auto partAssignment: partitionToProcessorAssignments) {
            for (auto partRange: partAssignment) {
                // First, collect the sizes of all partitions in partitionSizes
                // Only later compute the exclusive prefix sum.
                assert(partRange.part_id < partitionSizes.size());
                partitionSizes[partRange.part_id] += partRange.length;

                // Assume, that the partition ids are consecutive.
                _numPartitionsGlobal = std::max(_numPartitionsGlobal, partRange.part_id);

                // Also count the number of patterns to later infer the number of ReStore
                // blocks. Assume, that there are no overlapping partitions.
                numPatterns += partRange.length;
            }
        }
        _numPartitionsGlobal++;
        assert(_numPartitionsGlobal == _partitionIdOffset.size());
        assert(_numPartitionsGlobal == partitionSizes.size());

        // Compute the exclusive prefix sum over the partition sizes to obtain the
        // partition's id offset.
        assert(_partitionIdOffset.size() == partitionSizes.size());
        std::exclusive_scan(partitionSizes.begin(), partitionSizes.end(), _partitionIdOffset.begin(), 0);
        assert(_partitionIdOffset.size() == _numPartitionsGlobal);

        // Currently, the number of ReStore blocks is equal to the number of patterns across all partitions in the MSA.
        _numReStoreBlocks = numPatterns;

        const auto& localPartitionAssignment = partitionToProcessorAssignments[ParallelContext::rank_id()];
        _numPartitionsLocal                  = localPartitionAssignment.num_parts();
        assert(_numPartitionsLocal > 0);
        assert(_numPartitionsLocal <= _numPartitionsGlobal);
        return localPartitionAssignment;
    }
};