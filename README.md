# pyMorton
[![DOI](https://zenodo.org/badge/367796024.svg)](https://zenodo.org/badge/latestdoi/367796024)

pyMorton is a Python module with methods to efficiently map multidimensional data to a single dimension while preserving locality of the data points.

This mapping is commonly known as Z-order, Lebesgue curve, Morton space filling curve, Morton order or Morton code.

![](https://upload.wikimedia.org/wikipedia/commons/3/30/Z-curve.svg)
*Image by David Eppstein, 2008*

The Morton code of a multi-dimensional data point is calculated by bitwise interlacing the binary representations of its coordinate values.

pyMorton provides two functions for handling the encoding and decoding of data points with _arbitrary_ dimensionality and _arbitrary_ coordinate size:
 
```python
interlace(*data_point: int, dims: int = None, bits_per_dim: int = None) -> int
```
```python
deinterlace(code_point: int, dims: int = 3) -> List[int]
```

When handling large multi-dimensional dataset (n > 10.000), pyMorton offers some simple  but convenient means of parallelizing the Morton encoding and decoding:
 
```python
par_interlace(data_points: List[List[int]], dims: int = None, bits_per_dim: int = None) -> List[int]
```
```python
par_deinterlace(code_points: List[int], dims: int = 3) -> List[List[int]]
```

Given the Morton codes of a multi-dimensional dataset, we can perform multi-dimensional range search using only a one-dimensional data structure. 
For range searching, pyMorton offers two functions for calculating the necesaary `LITMAX` and `BIGMIN` values:
```python
prev_morton(code_point: int, rmin_code: int, rmax_code: int, dims: int = 3) -> int
```
```python 
next_morton(code_point: int, rmin_code: int, rmax_code: int, dims: int = 3) -> int
```

This implementation is based on the following paper 

> Tropf, Herbert, and Helmut Herzog. "Multidimensional Range Search in Dynamically Balanced Trees." ANGEWANDTE INFO. 2 (1981): 71-77.

and it makes heavy use of the excellent [gmpy2 module](https://gmpy2.readthedocs.io/en/latest/).

## Installation
```bash
pip install pyMorton
```

## Usage

### Basics 
````python
import pyMorton as pm 
````
imports the module.
```python
code = pm.interlace(2,16,8)
```
interlaces the 3D data point `(2,16,8)` into Morton code point `10248`.

When explicitly specify dimensionality and bits per dimension of your data point 
```python
code = pm.interlace(2,16,8, dims=3, bits_per_dim=5)
```
performance will benefit substantially. 

```python
pm.deinterlace(4711)
```
deinterlaces the Morton code point `4711` into the 3D data point `(29,1,3)`.

### Parallel interlacing/deinterlacing
Given a potentially large list of n-dimensional `data_points`
````python
from random import randrange

bit_size = 16
max_val = 2**bit_size - 1
no_samples = 10**6
data_points = [(randrange(0, max_val), randrange(0, max_val), randrange(0, max_val)) for i in range(no_samples)]
```` 
we can speed up things by using `par_interlace` and `par_deinterlace`
```python
morton_codes = pm.par_interlace(data_points, dims=3, bits_per_dim=16)
data_points == par_deinterlaces(morton_codes, dims=3)
````

### Range searching   
![](https://i.postimg.cc/qRQfSY80/tropf-figure-9.png)
*Image by Tropf and Herzog, 1981*


When range searching, we can prune the search space by calculating `BIGMIN` (aka "GetNextZ-address") and `LITMAX` (aka "GetPrevZ-address") values.     
```python
point = pm.interlace(6, 3, dims=2)  # => 30
rmin = pm.interlace(5, 3, dims=2)   # => 27
rmax = pm.interlace(10, 5, dims=2)  # => 102

BIGMIN = pm.next_morton(point, rmin, rmax, dims=2) # => 31
LITMAX = pm.prev_morton(point, rmin, rmax, dims=2) # => 27
```
In addition, we can easily check if a given Morton code point is within a specified range
```python 
pm.in_range(58,27,102, dims=2) # => False
pm.in_range(49,27,102, dims=2) # => True
```
 
## Citation
```bibtex
@misc{rmrschub_2021_pyMorton,
    author       = {René Schubotz},
    title        = {{pyMorton: Multi-dimensional indexing using Morton space filling curves.}},
    month        = may,
    year         = 2021,
    doi          = {10.5281/zenodo.4771019},
    version      = {0.0.1},
    publisher    = {Zenodo},
    url          = {https://github.com/rmrschub/pyMorton}
    }
```



## License
<a rel="license" href="http://creativecommons.org/licenses/by-nc-sa/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by-nc-sa/4.0/80x15.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by-nc-sa/4.0/">Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License</a>.