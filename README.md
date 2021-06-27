# ClusterCalibration
Cluster Calibration for the Atlas Detector

This project was meant to recalibrate the deposited energy in the atlas calorimeters.

## HPCC Network
To run this data 

```bash
# if partition is matador
sbatch matador_job.sh 
# if partition is quanah
sbatch quanah_job.sh 
```

## Local Machines
To run this data

```bash
./run.sh
```


When the code is finished running it will put all the results in a new folder.

If you have any questions don't hesitate to contact me at "madison.howard@ttu.edu"
