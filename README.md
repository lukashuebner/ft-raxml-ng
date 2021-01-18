# Fault-Tolerant RAxML Next Generation (FT-RAxML-NG)

## Introduction

RAxML-NG is a phylogenetic tree inference tool which uses maximum-likelihood (ML) optimality criterion. Its search heuristic is based on iteratively performing a series of Subtree Pruning and Regrafting (SPR) moves, which allows to quickly navigate to the best-known ML tree. RAxML-NG is a successor of RAxML (Stamatakis 2014) and leverages the highly optimized likelihood computation implemented in [*libpll*](https://github.com/xflouris/libpll) (Flouri et al. 2014).

RAxML-NG offers improvements in speed, flexibility and user-friendliness over the previous RAxML versions. It also implements some of the features previously available in ExaML (Kozlov et al. 2015), including checkpointing and efficient load balancing for partitioned alignments (Kobert et al. 2014).

## Fault Tolerance
FT-RAxML-NG is expanded to include support for failing nodes during computations on large clusters. [ULFM](https://fault-tolerance.org/) is used to detect failures and fix the MPI communicator. ULFM is a expansion of the OpenMPI MPI implementation. To get failure tolerance, you need to compile FT-RAxML-NG against ULFM and run in using ULFM's `mpirun`. FT-RAxML-NG will then automatically detect and mitigate failing ranks.

### Installing ULFM on a cluster
If your cluster environment does not include an up-to date ULFM installation, you can try installing one into your home directory.

0. First, create the target directry, for example `mkdir "$HOME/ulfm"`.
1. Download and extract the current versionof ULFM from [the official repo](https://bitbucket.org/icldistcomp/ulfm2/downloads/)
2. Configure ULFM using `./configure --prefix="$HOME/ulfm"`
3. Build using `make -j`
4. Install ULFM locally using `make install`

Ensure to set the following environment variables when _compiling_ and _running_ FT-RAxML-NG. 
```
# In this example $HOME/ulfm" is the prefix we installed ULFM into
export PATH="$HOME/ulfm/bin:$PATH"
export CPATH="$HOME/ulfm/src:$CPATH"
export LD_LIBRARY_PATH="$HOME/ulfm/lib:$LD_LIBRARY_PATH"
```

Your can verify that ULFM is used instead of a pre-installed MPI implementation by checking the output of `which mpirun`, which should point to your just installed binary. You can then proceed to build the MPI version of FT-RAxML-NG (see below).

## Installing (FT-)RAxML-NG

Please clone this repository and build RAxML-NG from scratch.

1. **Install the dependecies.** On Ubuntu (and other Debian-based systems), you can simply run:
```
sudo apt-get install flex bison libgmp3-dev
```
For other systems, please make sure you have following packages/libraries installed:  
[`GNU Bison`](http://www.gnu.org/software/bison/) [`Flex`](http://flex.sourceforge.net/) [`GMP`](https://gmplib.org/)

2. **Build FT-RAxML-NG.**

MPI version:

```
git clone --recursive https://github.com/amkozlov/raxml-ng
cd raxml-ng
mkdir build && cd build
cmake -DUSE_MPI=ON ..
make
```

## Running FT-RAxML-NG
Make sure you are using the mpirun shipped with ULFM. You can test this using `which mpirun`.  We recommend the following settings to reduce the number of false positive failure reports:
```
mpirun -n $SLURM_NTASKS \
  --mca mpi_ft_enable true \
  --mca mpi_ft_detector_thread true \
  --mca mpi_ft_detector_period 0.3 \
  --mca ft_detector_timeout 1 \
  ../raxml-ng-mpi --search --msa inputMSA.phy --model inputModel.model --threads 1
```
This will increase ULFM's heartbeat period to 300 ms, and timeout to 1 s. It will also enable a separate thread responsible for sending and receiving heartbeats. See [this](https://fault-tolerance.org/2020/01/21/spurious-errors-lack-of-mpi-progress-and-failure-detection/) article on the ULFM website for details.

### Limitations
Some RAxML-ng features are not failure-tolerant yet. For example, currently only `-search` mode without bootstrap replicas is supported. Also, only fine-grained parallelization with a single thread per rank is supported. 

### Numerical Instability of Allreduce Operations
Allreduce operations on floating-point values are numerically unstable. If the number of PEs
which take part in the allreduce operation changes, the result might change as well. This is,
because floating-point operations are only approximately associative and commutative. The
changed order of operations will cause the small inaccuracies to pile up differently.
This has impacts on the reproducibility of tree searches. When no failure occurred, we
can always reproduce a result by using the same number of nodes, cores per node, SIMD kernel, compiler, linker, OS version and the same random seed on the same system.
If a failure occurred, RAxML-ng will conduct different allreduce operations with a different
number of PEs. To reproduce this result, we would have to simulate a failure at the exact same
moment in the tree search. Implementing either a numerically stable allreduce operation or a
failure-log to enhance reproducibility is subject of future work.

## Documentation and Support

RAxML-NG, including documentation, a tutorial, and support can be found [here](https://github.com/amkozlov/raxml-ng). 

## License and citation

The code is currently licensed under the GNU Affero General Public License version 3.

When using the fault-tolerant FT-RAxML-NG, please cite [this paper](https://www.biorxiv.org/content/10.1101/2021.01.15.426773v1):

Lukas Hübner, Alexey M. Kozlov, Demian Hespe, Peter Sanders, Alexandros Stamatakis (2021)
**Exploring Parallel MPI Fault Tolerance Mechanisms for Phylogenetic Inference with RAxML-NG**
*Preprint at bioRxiv*
doi:[10.1101/2021.01.15.426773](https://doi.org/10.1101/2021.01.15.426773)

Additionally, please mention that FT-RAxML-NG is an extension of RAxML-NG and cite [this paper](https://doi.org/10.1093/bioinformatics/btz305) as well:

Alexey M. Kozlov, Diego Darriba, Tom&aacute;&scaron; Flouri, Benoit Morel, and Alexandros Stamatakis (2019)
**RAxML-NG: A fast, scalable, and user-friendly tool for maximum likelihood phylogenetic inference.** 
*Bioinformatics, btz305* 
doi:[10.1093/bioinformatics/btz305](https://doi.org/10.1093/bioinformatics/btz305)

## References

* Stamatakis A. (2014)
**RAxML version 8: a tool for phylogenetic analysis and post-analysis of large phylogenies.**
*Bioinformatics*, 30(9): 1312-1313.
doi:[10.1093/bioinformatics/btu033](http://dx.doi.org/10.1093/bioinformatics/btu033)

* Flouri T., Izquierdo-Carrasco F., Darriba D., Aberer AJ, Nguyen LT, Minh BQ, von Haeseler A., Stamatakis A. (2014)
**The Phylogenetic Likelihood Library.**
*Systematic Biology*, 64(2): 356-362.
doi:[10.1093/sysbio/syu084](http://dx.doi.org/10.1093/sysbio/syu084)

* Kozlov A.M., Aberer A.J., Stamatakis A. (2015)
**ExaML version 3: a tool for phylogenomic analyses on supercomputers.**
*Bioinformatics (2015) 31 (15): 2577-2579.*
doi:[10.1093/bioinformatics/btv184](https://doi.org/10.1093/bioinformatics/btv184)

* Kobert K., Flouri T., Aberer A., Stamatakis A. (2014)
**The divisible load balance problem and its application to phylogenetic inference.**
*Brown D., Morgenstern B., editors. (eds.) Algorithms in Bioinformatics, Vol. 8701 of Lecture Notes in Computer Science. Springer, Berlin, pp. 204–216*
