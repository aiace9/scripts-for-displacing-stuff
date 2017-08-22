OUT_FILE = """ &control
    calculation = '{calculation}',
    prefix = '{prefix}',
    restart_mode = 'from_scratch'
    pseudo_dir = '{pseudo_dir}',
    outdir = '{outdir}',
    wf_collect = .false.
    nstep = 250
    disk_io = 'none'
    verbosity  = 'high'
    etot_conv_thr = 1.0D-4
    forc_conv_thr = 5.0D-4
    max_seconds = {max_seconds}
    tstress = .true.
    tprnfor = .true.
 /
 &system
    input_dft = "vdW-DF"
    ibrav =  {crystal_structure},
    nat =  {nat},
    ntyp = 4,
    ecutwfc = {e_cutoff}
    ecutrho = {e_cutroh}
    spline_ps = .true.
 /
 &ELECTRONS
    electron_maxstep = 250
    conv_thr    = 1.D-8,
    mixing_beta = 0.5D0,
 /
 &IONS
   ion_dynamics='bfgs'
!  trust_radius_max = 0.2
 /
 &CELL
    cell_dynamics='bfgs'
 /

ATOMIC_SPECIES
{ATOMIC_SPECIES}

CELL_PARAMETERS {cell}
 {x}
 {y}
 {z}

ATOMIC_POSITIONS {{crystal}}
{ATOMIC_POSITIONS}
K_POINTS (automatic)
  {K_points}
"""

JOB_FILE = """#PBS -N {processName}
#PBS -q regular
#PBS -l nodes=1:ppn=20
#PBS -l walltime={wallTime}

NCPU=`wc -l < $PBS_NODEFILE`
NNODES=`uniq $PBS_NODEFILE | wc -l`

NOW=`date +%H:%M-%a-%d/%b/%Y`

echo ------------------------------------------------------
echo 'This job is allocated on '$NCPU' cpu(s)'
echo 'Job is running on node(s): '
cat  $PBS_NODEFILE
echo ------------------------------------------------------
echo WORKINFO:
echo PBS: the job was submitted     $NOW
echo PBS: qsub is running on        $PBS_O_HOST
echo PBS: originating queue is      $PBS_O_QUEUE
echo PBS: executing queue is        $PBS_QUEUE
echo PBS: working directory is      $PBS_O_WORKDIR
echo PBS: current home directory is $PBS_O_HOME
echo ""
echo JOBINFO:
echo PBS: execution mode is    $PBS_ENVIRONMENT
echo PBS: job identifier is    $PBS_JOBID
echo PBS: job name is          $PBS_JOBNAME
echo ""
echo NODEINFO:
echo PBS: node file is         $PBS_NODEFILE
echo PBS: number of nodes is   $NNODES
echo PBS: number of cpu        $NCPU
echo ------------------------------------------------------
cd $PBS_O_WORKDIR
module load openmpi/1.8.3/intel/14.0 mkl/11.1 intel/14.0
mpirun -np $NCPU ~/bin/pw.x -in {fFile} > {out}
"""
