How to run GenotypeRecurrentCNVs_tutorial.ipynb on Expanse

1. Copy GenotypeRecurrentCNVs_tutorial.ipynb to a working directory

2. Start an interactive job; make sure to give it at least 10G of memory 
   (I used 20G to be safe)
    srun --partition=ind-shared --pty --account=ddp195 --time 04:00:00 \
    --nodes=1 --mem 20G --cpus-per-task 1 --ntasks-per-node=1  --export=ALL /bin/bash

3. Create and activate a conda environment containing jupyterlab and pandas (you can use mine)
    conda env create -f jupyterEnv.yaml
    conda activate jupyterEnv
    -- OR -- use mine on Expanse
    conda activate /expanse/projects/sebat1/msacks/Utils/jupyterEnv/

4. Launch Jupyter Lab from this compute node (NOT your login node)
    jupyter notebook --no-browser --port=8888 --ip=0.0.0.0

5. On your laptop/ local machine, open a new terminal window

6. Make an SSH tunnel (I used port 8890 on my laptop) filling
   in the appropriate values and enter your credentials 
   If there's no output then the command is working.
    ssh -N -L 8890:<compute node>:8888 <username>@login.<your_cluster>.edu
   Here is the command I ran on my laptop based on the interactive job I started:
    ssh -N -L 8890:exp-15-01:8888 msacks@login.expanse.sdsc.edu

7. Go to your browser on your laptop/ local machine (I used Chrome) and go to:
    localhost:8890/tree

8. Navigate to your working directory and open up GenotypeRecurrentCNVs_tutorial.ipynb

You should now able to run the tutorial from your browser! 
This essentially runs the notebook in an interactive job then sets 
up a tunnel to your local machine, where you can run the notebook
in your browser.