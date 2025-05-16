# WPSD.jl

Code in Julia for computing 
[well-separated pairs decomposition](https://en.wikipedia.org/wiki/Well-separated_pair_decomposition), and
approximating the diameter of a point set.

Installation
------------
Uses some Julia code from the Frechet distance Julia
[Implementation](https://github.com/sarielhp/FrechetDist.jl). Currently
this is a bit of a hack - you need to create a symlink from
FrechetDist/src/ directory to *cg*.

For fun, I also implemented the diameter approximation algorithm from 
[here](https://sarielhp.org/p/00/diam.html) which is surprisingly easy
(and fast too?).

This is currently a bit of hack. I would release it as an official
package if there is interest.

Results
-------

| **Dim** | **N** | **RT Approx** | **RT Exact** | **Approx** |
|--------:|------:|--------------:|-------------:|-----------:|
| 3       | 2     | 0.000013      | 0.000000     | 1.00       |
| 3       | 4     | 0.000020      | 0.000000     | 1.00       |
| 3       | 8     | 0.000013      | 0.000000     | 1.00       |
| 3       | 16    | 0.000026      | 0.000000     | 0.98       |
| 3       | 32    | 0.000066      | 0.000001     | 1.00       |
