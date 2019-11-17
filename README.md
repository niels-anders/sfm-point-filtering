# Filtering of SfM point clouds

Code from the paper:
Anders N, Valente J, Masselink R, Keesstra S., 2019. *Comparing Filtering Techniques for Removing Vegetation from UAV-Based Photogrammetric Point Clouds*. Drones 3(61). [doi:10.3390/drones3030061](http://dx.doi.org/10.3390/drones3030061).

[Link to paper](https://www.mdpi.com/2504-446X/3/3/61)

## Installation

- dependencies: gdal-bin (tested with 2.2.3 from ubuntu repository)
- create a new python3 virtual environment with e.g.:

```bash
python3 -m venv venv
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

- folder='location to folder with las files'. Default=`example-data/`
- res=resolution of output dtm in meters. Default=`0.25`
- vi=vegetation index ('vvi' and 'exg' are supported, i.e. visible vegetation index and excessive greeness respectively). Default=`exg`
- threshold=vi threshold, default=`5e-5` for exg and `0.01` for vvi
- plot=True/False. Default=`False`
