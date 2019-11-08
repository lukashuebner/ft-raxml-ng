#!/usr/bin/env python3

import os

PE_PER_NODE = 20

def writeSbatch(jobDir, sbatchFile, numNodes, ranksPerNode, threadsPerRank, inputFile):
    with open(sbatchFile, "w") as f:
        f.write("#!/bin/bash\n")
        f.write("#SBATCH --nodes=%d\n" % numNodes)
        f.write("#SBATCH --ntasks-per-node=%d\n" % ranksPerNode)
        f.write("#SBATCH -t 01:00:00\n")
        f.write("#SBATCH -p normal\n")
        f.write("#SBATCH --mem=63720mb\n")
        f.write("#SBATCH --chdir=%s\n" % jobDir)
        f.write("#SBATCH -o out.txt\n")
        f.write("#SBATCH -e err.txt\n")
        f.write("\n")
        f.write("mpirun -n " + str(numNodes * ranksPerNode) + " ../raxml-ng-mpi --search --msa ../" + inputFile +
                " --model GTR+R --tree pars{10} --threads " + str(threadsPerRank) + " --redo\n")

def scheduleJob(sbatchFile):
    os.system("sbatch " + sbatchFile)

def run(dataset, numNodes, ranksPerNode, threadsPerRank):
    jobDir = "%s_%d@%d" % (dataset, ranksPerNode, numNodes)
    if (not os.path.isdir(jobDir)):
        os.mkdir(jobDir)
    sbatchFile = "%s.sbatch" % jobDir
    writeSbatch(jobDir, sbatchFile, numNodes, ranksPerNode, threadsPerRank, dataset + ".phy")
    #scheduleJob(sbatchFile)

def main():
    # DNA vs Amino Acids
    run(dataset = "dna_rokasD4", numNodes = 8, ranksPerNode = 20, threadsPerRank = 1)
    run(dataset = "aa_rokasA4", numNodes = 8, ranksPerNode = 20, threadsPerRank = 1)

    # Single part <-> Multi part
    run(dataset = "dna_rokasD1", numNodes = 18, ranksPerNode = 20, threadsPerRank = 1)
    run(dataset = "dna_rokasD2b", numNodes = 18, ranksPerNode = 20, threadsPerRank = 1)

    # Few patterns <-> many patterns    
    run(dataset = "dna_ShiD9", numNodes = 1, ranksPerNode = 20, threadsPerRank = 1)
    run(dataset = "dna_rokasD2a", numNodes = 1, ranksPerNode = 20, threadsPerRank = 1)
    
    # New nodes <-> many nodes
    run(dataset = "dna_rokasD1", numNodes = 4, ranksPerNode = 20, threadsPerRank = 1)
    run(dataset = "dna_rokasD1", numNodes = 30, ranksPerNode = 20, threadsPerRank = 1)
    
    # using pthreads <-> not using pthreads
    run(dataset = "dna_rokasD6", numNodes = 10, ranksPerNode = 20, threadsPerRank = 1)
    run(dataset = "dna_rokasD6", numNodes = 10, ranksPerNode = 1, threadsPerRank = 10)


if __name__ == "__main__":
    main()