# How to Run GenotypeRecurrentCNVs_tutorial.ipynb on an HPC Cluster

## 1. Copy the notebook

Copy `GenotypeRecurrentCNVs_tutorial.ipynb` to your working directory.

## 2. Start an interactive job

Make sure to allocate at least 10G of memory (20G recommended):

```bash
srun --partition=ind-shared --pty --account=<your_account> --time 04:00:00 \
    --nodes=1 --mem 20G --cpus-per-task 1 --ntasks-per-node=1 --export=ALL /bin/bash
```

## 3. Create and activate a conda environment

The environment needs `jupyterlab` and `pandas`:

```bash
conda env create -f jupyterEnv.yaml
conda activate jupyterEnv
```

Or activate an existing environment if available:

```bash
conda activate /path/to/your/jupyterEnv/
```

## 4. Launch Jupyter Lab

From the **compute node** (not your login node):

```bash
jupyter notebook --no-browser --port=8888 --ip=0.0.0.0
```

## 5. Set up SSH tunnel

On your local machine, open a new terminal window and create an SSH tunnel:

```bash
ssh -N -L <local_port>:<compute_node>:8888 <username>@login.<your_cluster>.edu
```

For example:

```bash
ssh -N -L 8890:exp-15-01:8888 username@login.cluster.edu
```

If there's no output, the command is working.

## 6. Access the notebook

1. Open your browser and go to: `localhost:<local_port>/tree`
2. Navigate to your working directory
3. Open `GenotypeRecurrentCNVs_tutorial.ipynb`

You should now be able to run the tutorial from your browser! This setup runs the notebook in an interactive job and tunnels it to your local machine.
