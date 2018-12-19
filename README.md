# blctn-input
Generate input files for Blue Action project.

```
➜ blctn-input.py --help
usage: blctn-input.py [-h] [-s STARTYEAR] [-e ENDYEAR] [-f]

optional arguments:
  -h, --help            show this help message and exit
  -s STARTYEAR, --startyear STARTYEAR
                        startyear to use for processing files [default 1979]
  -e ENDYEAR, --endyear ENDYEAR
                        stop processing in the beginning of endyear [default
                        2016]
  -f, --force           force overwriting existing files
```

## experiment 1
Perform SST and SIC adjustments as per instructions:

    1. Set minimum SST to -1.8 degC
    2. Set SST to -1.8 degC if SIC>0.9
    3. If SST>5 degC, set SIC to 0
    4. If SIC<0.9, we calculate SSTmax, where SSTmax=9.328*(0.729-SIC^3)-1.8. If SST>SSTmax, reduce SIC so that SST=SSTmax.

## experiment 2
Same as experiment 1, but with daily sea ice climatology. Redo the SST and SIC adjustments.

Perform SST and SIC adjustments as per instructions:

    1. Set minimum SST to -1.8 degC
    2. Set SST to -1.8 degC if SIC>0.9
    3. If SST>5 degC, set SIC to 0
    4. If SIC<0.9, we calculate SSTmax, where SSTmax=9.328*(0.729-SIC^3)-1.8. If SST>SSTmax, reduce SIC so that SST=SSTmax.

## experiment 3
Same as experiment 1, but with PDO removed from SST. Redo the SST and SIC adjustments.

Perform SST and SIC adjustments as per instructions:

    1. Set minimum SST to -1.8 degC
    2. Set SST to -1.8 degC if SIC>0.9
    3. If SST>5 degC, set SIC to 0
    4. If SIC<0.9, we calculate SSTmax, where SSTmax=9.328*(0.729-SIC^3)-1.8. If SST>SSTmax, reduce SIC so that SST=SSTmax.

## experiment 4
Same as experiment 1, but with AMO removed from SST. Redo the SST and SIC adjustments.

Perform SST and SIC adjustments as per instructions:

    1. Set minimum SST to -1.8 degC
    2. Set SST to -1.8 degC if SIC>0.9
    3. If SST>5 degC, set SIC to 0
    4. If SIC<0.9, we calculate SSTmax, where SSTmax=9.328*(0.729-SIC^3)-1.8. If SST>SSTmax, reduce SIC so that SST=SSTmax.
