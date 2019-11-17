# Filtering of SfM point clouds

code from the paper: https://www.mdpi.com/2504-446X/3/3/61

## Installation

- dependencies: gdal-bin (tested with 2.2.3 from ubuntu repository)
- create a new python3 virtual environment with e.g.:

```bash
python3 -m virtualenv venv
source venv/bin/activate
```

- install requirements.txt with 

```bash
pip install -r requirements.txt
```

## Run

run example with:

```bash
python main.py filter=isl_vi
```
Optional filters are `isl` (iterative surface lowering), `vi` (vegetation index), and a combination of both `isl_vi`
Optional arguments:
- folder='location to folder with las files'
- res=resolution of output dtm
- vi=vegetation index ('vvi' and 'exg' are supported, i.e. visible vegetation index and excessive greeness respectively)
- threshold=vi threshold
- plot=True/False 
