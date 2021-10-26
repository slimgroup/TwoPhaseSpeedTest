# TwoPhaseSpeedTest

3D Example of 2 phase flow solvers based on OPM

Simulations based on 2 homogeneous permeability/porosity models (discretized in `16^3` grid and `256^3` grid) are in `CO2SIM16.DATA` and `CO2SIM.DATA`.

TO DO:

```bash
flow CO2SIM.DATA/CO2SIM16.DATA
```

A simulation based on a `256^3` cube section of Compass model is also provided. First, run the scripts in `TakeCompassSubsection` to generate the permeability from velocity. For simplicity, you could download the input file to OPM solvers by `wget`, provided in `DownloadFile.sh`.

After `PERMX.txt`, `PERMY.txt`, `PERMZ.txt` and `PORO.txt` are downloaded, run

```julia
flow CO2SIM256Compass.DATA
```