# SETIDatFiles.jl

A Julia package for reading and writing SETI DAT files as produced by
[turboSETI], [seticore], and [bliss].

[turboSETI]: https://github.com/UCBerkeleySETI/turbo_seti
[seticore]: https://github.com/lacker/seticore
[bliss]: https://github.com/UCBerkeleySETI/bliss

## SETI DAT files

SETI DAT files describe Doppler-drifting signals (*hits*) that were detected in
radio telescope observations.  They consist of a header section containing
metadata about the observation followed by a white-space delimited table
describing the parameters of the detected Doppler-drifting signals.  The format
originated in `turboSETI` and was carried on by successor programs `seticore`
and `bliss`.  SETI DAT files usually have a `.dat` extension.  Due to the
organic evolution of the SETI DAT file format, there is no formal
specification.  `SETIDatFiles.jl` attempts to be able to read as many SETI DAT
file variations as possible while writing SETI DAT files in a very consistent
form.

## Installation

You can use Julia's `Pkg` package to add `SETIDatFiles` to your currently active
project:

```jl
import Pkg
Pkg.add(url="https://github.com/david-macmahon/SETIDatFiles.jl")
```

## Types

`SETIDatFiles` and provides two types to describe the data in a SETI DAT file:
- `DatFileHeader` contains the metadata of the DAT file header
- `DatFileHit` contains the parameters of a single Doppler-drifting signal

### `DatFileHeader`

`DatFileHeader` has fields for the SETI DAT file header metadata.  Where
possible, field names mirror SIGPROC Filterbank header names.

| Field         | Type    | Description                                          |
|:--------------|:--------|:-----------------------------------------------------|
| `fbh5name`    | String  | Name of the file containing the observational data   |
| `source_name` | String  | Name of the file containing the observational data   |
| `tstart`      | Float64 | The start time of the observation (MJD)              |
| `src_raj`     | Float64 | Source right ascension (J2000, hours)                |
| `src_dej`     | Float64 | Source declination (J2000, degrees)                  |
| `tsamp`       | Float64 | Integration time of each time sample (seconds)       |
| `foff`        | Float64 | Frequency resolution in (MHz, but in file it is Hz!) |
| `maxdr`       | Float64 | Maximum drift rate searched (Hz/s)                   |
| `nsamps`      | Int     | Number of samples (aka time steps) in FBH5 file      |

### `DatFileHit`

`DatFileHit` has fields for the parameters of a single hit from a DAT file.
Field names mirror those given in the static text of the DAT file header.
`Index` and `Coarse_Channel_Number` are 1-based, but they are stored in file as
0-based.

| Field                   | Type    | Description                              |
|:------------------------|:--------|:-----------------------------------------|
| `Top_Hit_Number`        | Int     | Hit index                                |
| `Drift_Rate`            | Float64 | Drift rate (Hz/s)                        |
| `SNR`                   | Float64 | Signal-to-noise ratio                    |
| `Uncorrected_Frequency` | Float64 | (note 1)                                 |
| `Corrected_Frequency`   | Float64 | (note 1)                                 |
| `Index`                 | Int     | Fine channel of `freq_start` (note 2)    |
| `freq_start`            | Float64 | Frequency of hit at start of observation |
| `freq_end`              | Float64 | Frequency of hit at end of observation   |
| `SEFD`                  | Float64 | System equivalent flux density (note 3)  |
| `SEFD_freq`             | Float64 | SEFD frequency (note 3)                  |
| `Coarse_Channel_Number` | Int     | Coarse channel containing hit (note 2)   |
| `Full_number_of_hits`   | Int     | (note 4)                                 |

- Note 1: `Uncorrected_Frequency` and `Corrected_Frequency` are under-defined.
  Different tools output different information here.  Most store the same value
  as `freq_start` here, but `bliss` stores the midpoint frequency between
  `freq_start` and `freq_end`.

* Note 2: `Index` and `Coarse_Channel_Number` refer to coarse channel and fine
  channel within the coarse channel that correspond to `freq_start`.  These
  values are stored in the DAT file as 0-based index values, but in `DatFileHit`
  they are 1-based to match Julia's indexing convention.  If the number of fine
  channels per coarse channel is not known, a valid approach is to set
  `Coarse_Channel_Number` 1 and `Index` to the hit's fine channel within the
  file.

* Note 3: `SEFD` and `SEFD_freq` seem intended to store sensitivity data about
  the telescope during the observation, but this data is often not present at
  signal searching time so these fields are often unused and set to 0.0.

* Note 4: `Full_number_of_hits` is another under-defined field.  Some uses are
  the number of *proto-hits* clustered near this hit in the detection plane or
  the number of fine channels that the detected signal occupies (i.e.
  proportional to the bandwidth of the signal).

## Functions

`SETIDatFiles` provides the `readdat` and `writedat` functions to read and write
SETI DAT files.  `readdat` tries to be as lenient as possible so it can parse
as many SETI DAT file variations as possible.  `writedat` is very consistent in
what it outputs.  This can be used to canonicalize a DAT file found in the wild:

```jl
writedat("canonicalized.dat", readdat("wild.dat"; fix_freq_end=true))
```

---

```jl
readdat(src; fix_freq_end=false) -> (hdr, hits)
readdat(f, src; fix_freq_end=false) -> (hdr, f(hits))
```

Read a SETI DAT file from `src` (an `IO` or a filename) returning a
`DatFileHeader` and a `Vector` of `DatFileHit` objects.  If `fix_freq_end=true`
is passed, the `freq_end` values from the file will be ignored and instead be
calculated as:

```jl
freq_end = freq_start + Drift_Rate/1e6 * hdr.tsamp * hdr.nsamps
```

A function or Type may be passed as the first parameter to modify the hits
before returning.  This can be used, for example, to have the hits returned as
a `StructArray` or `DataFrame` (which are not dependencies of this package).
For example:

```jl
julia> using SETIDatFiles, StructArrays

julia> hdr, hits = readdat(StructArray, "voyager.dat")

julia> hits.Drift_Rate
3-element Vector{Float64}:
 -0.367353
 -0.367353
 -0.367353
 ```

---

```jl
writedat(dest, hdr, hits) -> nothing
writedat(dest, hdrhits)   -> nothing
```

Write a SETI DAT file to `dest` (an `IO` or a filename).  `hdr` is a
`DatFileHeader` and `hits` may be any `AbstractVector{DatFileHit}` such as
`Vector{DatFileHit}` or `StructVector{DatFileHit}`.  `hdr` and `hits` may be
given a separate arguments or may be passed as a tuple.  `writedat` returns
`nothing`.
