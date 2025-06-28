# CaloXDataAnalysis

Download the data file from [here](https://yofeng.web.cern.ch/yofeng/CaloX/)

Change the path in `data/datafiles.json`

First compile the functions in `utils/functions.cc`
```
root [0] .L utils/functions.cc+
```

DQM plots:

```
python3 prepareDQMPlots.py
python3 makeDQMPlots.py
```

Hodoscope for tracking
```
python3 makeHodoPlots.py
```

FERS vs DRS correlations
```
python3 checkFERSDRS.py
```
