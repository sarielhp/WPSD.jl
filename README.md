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

┌───────────┬────────┬────────────────┬───────────────┬──────────────┐
│ dimension │      N │ approx_runtime │ exact_runtime │ approx_ratio │
│     Int64 │  Int64 │         String │        String │       String │
├───────────┼────────┼────────────────┼───────────────┼──────────────┤
│         3 │      2 │       0.000013 │      0.000000 │         1.00 │
│         3 │      4 │       0.000024 │      0.000000 │         1.00 │
│         3 │      8 │       0.000021 │      0.000000 │         1.00 │
│         3 │     16 │       0.000026 │      0.000001 │         1.00 │
│         3 │     32 │       0.000058 │      0.000001 │         1.00 │
│         3 │     64 │       0.000142 │      0.000004 │         1.00 │
│         3 │    128 │       0.000346 │      0.000017 │         1.00 │
│         3 │    256 │       0.000589 │      0.000067 │         1.00 │
│         3 │    512 │       0.001254 │      0.000270 │         1.00 │
│         3 │   1024 │       0.002125 │      0.001080 │         1.00 │
│         3 │   2048 │       0.003475 │      0.004333 │         1.00 │
│         3 │   4096 │       0.005668 │      0.017476 │         1.00 │
│         3 │   8192 │       0.007498 │      0.068837 │         1.00 │
│         3 │  16384 │       0.010921 │      0.275865 │         1.00 │
│         3 │  32768 │       0.014172 │      1.108290 │         1.00 │
│         3 │  65536 │       0.022967 │      4.478142 │         1.00 │
│         3 │ 131072 │       0.036131 │     17.887796 │         1.00 │
└───────────┴────────┴────────────────┴───────────────┴──────────────┘

