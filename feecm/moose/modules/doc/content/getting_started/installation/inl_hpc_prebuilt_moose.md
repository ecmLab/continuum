# INL HPC Pre-Built MOOSE

While operating on one of the [!ac](INL) [!ac](HPC) clusters, there exists the option of using
pre-built versions of MOOSE. To request access to these clusters, please follow the instructions on
[INL's Nuclear Computational Resource Center](https://inl.gov/ncrc/) website.

Once access has been granted, log into Sawtooth or Lemhi using either [inl/hpc_ondemand.md]
Interactive Shell services, or directly by following our [SSH Primer](inl/hpc_remote.md).

!include installation/clone_moose.md

## Load Modules

Load the following modules:

```bash
module load use.moose moose-apps moose
```

!alert warning
If you receive an error about modules not being known, please make sure you are logged into either
Sawtooth or Lemhi.

Once loaded, `moose-opt` becomes available. You need now only provide input files to run
simulations. Example input files are available within the MOOSE repository (next section).

## Run Example

To run `moose-opt` using an example input file from within the MOOSE repository, perform the
following:

```bash
cd ~/projects/moose/examples/ex01_inputfile
moose-opt -i ex01.i
```

`moose-opt` should run to completion without error. A resulting exodus file is generated in the same
directory:

```pre
ex01_out.e
```

## View Results

You can use HPC OnDemand to view the results file remotely. Head on over to
[HPC OnDemand Dashboard](https://hpcondemand.inl.gov/pun/sys/dashboard), and select:
`Interactive Apps` and then `Linux Desktop with Visualization`. Next, select your cluster (such as
Sawtooth), the amount of time you believe you need, and then click `Launch`.

It may take some time before your 'Visualization Job' becomes available. When it does, simply click
on it, and you will be presented a [!ac](GUI) desktop within your web browser. From here, you can
open visualization applications (such as Paraview), and open your results file.

To use Paraview, open a terminal by clicking `Applications` at the top left, then click
`Terminal Emulator`. A terminal window will open. Enter the following commands:

```bash
module load paraview
paraview
```

Paraview should open. From here, you can select `File`, `Open`, and navigate to

```pre
~/projects/moose/examples/ex01_inputfile
```

You should see your results file listed (`ex01_out.e`). Double click this file and enjoy the
results!

In many cases it is more desirable to view results using your local machine. This is done by copying
result files from the remote machine to your local machine using `scp` or `rsnyc`.

!alert note title=Copying files from remote HPC machine to your machine
Copying files from an HPC cluster to your machine first requires that you follow the instructions:
[inl/hpc_remote.md]. This workflow is for advanced users, comfortable with modifying their SSH
settings, and familiar with their terminal in general.

Example commands you would use with a terminal to copy files from HPC to your local machine's
Downloads folder:

```bash
cd ~/Downloads
scp <your hpc user id>@hpclogin:~/projects/moose/examples/ex01_inputfile/ex01_out.e .
```
