# Snippets GMT6 (modified) 

Snippets GMT6 provides GMT6 scripts to visualize [FFI](https://github.com/SELT-tsukuba/FFI) results (e.g. fort.40) obtained under `ICMN=1,2,5`.


## Lineups

- [Statics](#statics)
  - [ffi_solution.bash](#ffi_solutionbash)
  - [ffi_MapMeca.bash](#ffi_MapMecabash)
  - [ffi_MapX.bash](#ffi_MapXbash)
  - [ffi_RupStrDip.bash](#ffi_RupStrDipbash)
  - [ffi_NDC.bash](#ffi_NDCbash)

- [Dynamics](#dynamics)
  - [ffi_PlaneMecaSnap.bash](#ffi_PlaneMecaSnapbash)
  - [ffi_MapMecaSnap.bash](#ffi_MapMecaSnapbash)
  - [ffi_MapXSnap.bash](#ffi_MapXSnapbash)

- [Waveform fittings](#waveform-fittings)
  - [ffi_WaveComp.bash](#ffi_wavecompbash)
  - [ffi_AziEquiStaMap.bash](#ffi_aziequistamapbash)


## Installing

1. Open bash on your computer and navigate to the directory where you want to put this repository (Your computer should already be setup with Git and a bash shell interface. If not, please refer to the [Guidance](https://github.com/SELT-tsukuba/Guidance) before continuing).

2. The `git clone` command copies this repository from GitHub to your local machine.
```
$ git clone https://github.com/SELT-tsukuba/SnippetsGMT6_m.git
```

3. Once you have this repository, go to `SnippetsGMT6_m/src` and try `make`.

4. Finally, add a `PATH` to `SnippetsGMT6_m/bin` in `~/.bashrc` of your local machine.
```
export PATH=~/ Path To /SnippetsGMT6_m/bin:$PATH
```


## Description


### Statics


#### ffi\_solution.bash
<img src="./images/ffi_solution.png" width="500px">

```
$ ffi_solution.bash
$ open ffi_solution.pdf
```


#### ffi\_MapMeca.bash
<img src="./images/ffi_MapMeca.png" width="500px">

```
$ ffi_MapMeca.bash
$ open ffi_MapMeca.png
```

#### ffi\_MapX.bash
<img src="./images/ffi_MapX.png" width="500px">

```
$ ffi_MapX.bash  [ s p t b ] # s: strike / p: P axis / t: T axis / b: B axis
$ open ffi_MapX.png
```


#### ffi_RupStrDip.bash
<img src="./images/ffi_RupStrDip.png" width="500px">

```
$ ffi_RupStrDip.bash  [total duration]
$ open ffi_RupStrDip.png
```


#### ffi\_NDC.bash
<img src="./images/ffi_NDC.png" width="500px">

```
$ ffi_NDC.bash [snap interval] [total duration] [Average: 0,  Moment: 1]
$ open ffi_NDC.png
```

#### MapModelPlane.bash
<!--img src="./images/MapModelPlane.png" width="500px"-->

```
$ MapModelPlane.bash [strike] [dip] [lat] [lon] [depth] [xx] [yy] [mn] [nn] [m0] [n0]
$ open MapModelPlane.png
```

### Dynamics


#### ffi\_PlaneMecaSnap.bash
<img src="./images/ffi_PlaneMecaSnap.png" width="500px">

```
$ ffi_PlaneMecaSnap.bash  [snap interval] [total duration] [Average: 0,  Moment: 1] [Rotate-Focal-Mechanism : 0, No-Rotation : 1]
$ open ffi_PlaneMecaSnap.png
```


#### ffi_MapMecaSnap.bash
<img src="./images/ffi_MapMecaSnap.png" width="500px">

```
$ ffi_MapMecaSnap.bash [snap interval] [total duration] [Average: 0,  Moment: 1]
$ open ffi_MapMecaSnap.pdf
```

#### ffi\_MapXSnap.bash 
<img src="./images/ffi_MapXSnap.png" width="500px">

```
$ ffi_MapXSnap.bash  [snap interval] [total duration] [Average: 0,  Moment: 1] [ s p t b ] 
```


### Waveform fittings


#### ffi\_WaveComp.bash
<img src="./images/ffi_WaveComp.png" width="500px">

```
$ ffi_WaveComp.bash
$ open ffi_WaveComp.pdf
```

#### ffi\_AziEquiStaMap.bash
<img src="./images/ffi_AziEquiStaMap.png" width="500px">

```
$ ffi_AziEquiStaMap.bash
$ open ffi_AziEquiStaMap.pdf
```


## Plot all

```
$ ffi_plot_all.bash
$ open ffi_*.png
```


## File contents

### Statics

- `pdtdis.dat`
  ```
  1  2  3   4   5    6    7      8     9    10   11   12   13   14   15     16    17    18    19    20     21     22      23      24      25      26      27      28
  n, m, dx, dy, lat, lon, depth, slip, Mrr, Mss, Mee, Mrs, Mre, Mse, nexpo, str1, str2, dip1, dip2, rake1, rake2, trendp, trendt, trendb, plungp, plungt, plungb, NDC
  ```

- `pddis.dat`
  ```
  1    2    3   4   5     6
  lat, lon, dx, dy, slip, depth
  ```

- `faultline.dat`
  ```
  lat, lon, depth
  ```

- `tpdt.dat`
  ```
  Mrr, Mss, Mee, Mrs, Mre, Mse
  str1, dip1, rake1, str2, dip2, rake2
  sm, NDC, 0., 0., 0., 0.
  ```

- `mrf.dat`
  ```
  1     2           3    4    5    6    7    8
  time, momentrate, Mrr, Mss, Mee, Mrs, Mre, Mse
  ```

- `knotplot.dat`
  ```
  1  2  3    4    5      6       7    8     9   10
  n, m, lat, lon, depth, strike, dip, rake, dx, dy
  ```

- `knotcorner.dat`
  ```
  1  2  3       4       5         6       7       8         9       10      11        12      13      14        15      16   17    18  19
  n, m, lat_bl, lon_bl, depth_bl, lat_tl, lon_tl, depth_tl, lat_tr, lon_tr, depth_tr, lat_br, lon_br, depth_br, strike, dip, rake, dx, dy
  ```

- `knotedge.dat`
  ```
  lat, lon, depth
  ```


### Dynamics

- `snap.dat`
  ```
  1  2  3   4   5   6         7    8    9    10   11   12   13    14    15    16    17     18     19   20   21     22      23      24      25      26      27      28
  n, m, tw, dx, dy, sliprate, Mrr, Mss, Mee, Mrs, Mre, Mse, str1, str2, dip1, dip2, rake1, rake2, lat, lon, depth, trendp, trendt, trendb, plungp, plungt, plungb, NDC
  ```

- `snap2.dat`
  ```
  1  2  3   4   5   6         7    8    9
  n, m, tw, dx, dy, sliprate, lat, lon, depth
  ```

- `tw_mec_xy.dat`
  ```
  1   2   3   4    5    6    7    8    9    10    11    12    13    14     15     16   17   18
  tw, dx, dy, Mrr, Mss, Mee, Mrs, Mre, Mse, str1, str2, dip1, dip2, rake1, rake2, lat, lon, depth
  ```

- `faultline.dat`
  ```
  lat, lon, depth
  ```

- `vr_str.dat` or `vr_dip.dat`
  ```
  location, time, sliprate
  ```


## Acknowledgments

- The [GNU Fortran](https://gcc.gnu.org/fortran/) portion of the Snippets GMT6 is based on the [PyD](https://github.com/SELT-tsukuba/PyD).
# SnippetsGMT6_m
