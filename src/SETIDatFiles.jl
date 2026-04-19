module SETIDatFiles

export DatFileHeader, DatFileHit, readdat, writedat

using DelimitedFiles, Printf, AstroAngles, PrecompileTools

const SEPARATOR = "# -------------------------- o --------------------------"

"""
    DatFileHeader

Structure for SETI DAT file header metadata.  Where possible, field names mirror
SIGPROC Filterbank header names.

| Field      | Type    | Description                                          |
|:-----------|:--------|:-----------------------------------------------------|
| `fbh5name` | String  | Name of the file containing the observational data   |
| `tstart`   | Float64 | The start time of the observation (MJD)              |
| `src_raj`  | Float64 | Source right ascension (J2000, hours)                |
| `src_dej`  | Float64 | Source declination (J2000, degrees)                  |
| `tsamp`    | Float64 | Integration time of each time sample (seconds)       |
| `foff`     | Float64 | Frequency resolution in (MHz, but in file it is Hz!) |
| `maxdr`    | Float64 | Maximum drift rate searched (Hz/s)                   |
| `nsamps`   | Int     | Number of samples (aka time steps) in FBH5 file      |
"""
struct DatFileHeader
    fbh5name::String
    source_name::String
    tstart::Float64
    src_raj::Float64
    src_dej::Float64
    tsamp::Float64
    foff::Float64
    maxdr::Float64
    nsamps::Int
end

"""
    DatFileHit

Structure for a hit from a DAT file.  Field names mirror DAT file names.
`Index` and `Coarse_Channel_Number` are 1-based, but in file they are 0-based.

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

The extended help contains the notes.

# Extended help

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
"""
struct DatFileHit
    #=  1 =# Top_Hit_Number::Int
    #=  2 =# Drift_Rate::Float64
    #=  3 =# SNR::Float64
    #=  4 =# Uncorrected_Frequency::Float64
    #=  5 =# Corrected_Frequency::Float64
    #=  6 =# Index::Int
    #=  7 =# freq_start::Float64
    #=  8 =# freq_end::Float64
    #=  9 =# SEFD::Float64
    #= 10 =# SEFD_freq::Float64
    #= 11 =# Coarse_Channel_Number::Int
    #= 12 =# Full_number_of_hits::Int
end

DatFileHit(v::AbstractVector) = DatFileHit(v...)

"""
    readdat(src; fix_freq_end=false) -> (hdr, hits)
    readdat(f, src; fix_freq_end=false) -> (hdr, f(hits))

Read a SETI DAT file from `src` (an `IO` or a filename) returning a
`DatFileHeader` and a `Vector` of `DatFileHit` objects.  If `fix_freq_end=true`
is passed, the `freq_end` values from the files will be ignored and instead be
calculated as:

```jl
freq_end = freq_start + Drift_Rate * hdr.tsamp * hdr.nsamps
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
"""
function readdat(io::IO; fix_freq_end=false)
    # Read header
    readline(io) # separator
    fnameline = readline(io)
    readline(io) # separator
    srcline = readline(io) # src
    mjdradecline = readline(io) # mjd/ra/dec
    paramsline = readline(io) # params
    readline(io) # separator
    readline(io) # field names
    readline(io) # separator

    # Parse fname line
    fnameline = replace(fnameline, ":"=>" "; count=1)
    fnameline = replace(fnameline, r"\s+$"=>"") # trim trailing whitespace
    fbh5name = last(split(fnameline, r"\s+"; limit=4))
    
    # Parse source line
    srcline = replace(srcline, ":"=>" "; count=1)
    src = last(split(srcline, r"\s+"; limit=3))

    # Parse MJD/RA/Dec line
    mjdradecline = replace(mjdradecline, r"(MJD|RA|DEC):"=>s"\1 ")
    words = split(mjdradecline, r"\s+")
    mjd = length(words) < 3 ? 51544.0 :
        something(tryparse(Float64, words[3]), 51544.0)
    ra = length(words) < 5 ? 0.0 :
        try
            hms2ha(words[5])
        catch
            0.0
        end
    dec = length(words) < 7 ? 0.0 :
        try
            dms2deg(words[7])
        catch
            0.0
        end

    # Parse DELTAF/DELTAT/MAX_DR/OBSLENGTH line
    paramsline = replace(paramsline, ":"=>" ")
    words = split(paramsline, r"\s+")
    tsamp  = length(words) < 3 ? 1.0 : something(tryparse(Float64, words[3]), 1.0)
    foff   = length(words) < 5 ? 1.0 : something(tryparse(Float64, words[5]), 1.0)/1e6 # Store as MHz
    maxdr  = length(words) < 7 ? 0.0 : something(tryparse(Float64, words[7]), 0.0)
    obslen = length(words) < 9 ? 1.0 : something(tryparse(Float64, words[9]), 1.0)
    nsamps = round(Int, obslen/tsamp)

    hdr = DatFileHeader(
        fbh5name,
        src,
        mjd,
        ra,
        dec,
        tsamp,
        foff,
        maxdr,
        nsamps
    )

    # Read data (if any)
    hits = if eof(io)
        DatFileHit[]
    else
        dat = readdlm(io; comments=false)

        # Convert indices to 1-based
        dat[:,  6] .+= 1 # Index
        dat[:, 11] .+= 1 # Coarse_Channel_Number

        # Some DAT file writers output freq_end values that are not the ending
        # frequeuncy of the drift line.  For example, seticore outputs the same
        # value for freq_start and freq_end regardless of drift rate.  Passing
        # `fix_freq_end=true` will recompute the freq_end values based on
        # freq_start, Drift_Rate, and observation length.
        if fix_freq_end
            obslength = hdr.tsamp * hdr.nsamps
            # freq_end = freq_start + Drift_Rate * obslength
            dat[:, 8] .= dat[:, 7] .+ dat[:, 2] .* obslength
        end

        DatFileHit.(eachrow(dat))
    end

    (; hdr, hits)
end

function readdat(datname::AbstractString; fix_freq_end=false)
    #(; datname, open(readdat, datname)...)
    open(io->readdat(io; fix_freq_end), datname)
end

function readdat(f, src; fix_freq_end=false)
    h, v = readdat(src; fix_freq_end)
    h, f(v)
end

"""
    writedat(dest, hdr, hits) -> nothing
    writedat(dest, hdrhits)   -> nothing

Write a SETI DAT file to `dest` (an `IO` or a filename).  `hdr` is a
`DatFileHeader` and `hits` may be any `AbstractVector{DatFileHit}` such as
`Vector{DatFileHit}` or `StructVector{DatFileHit}`.  `hdr` and `hits` may be
given a separate arguments or may be passed as a tuple.  `writedat` returns
`nothing`.
"""
function writedat(io, hdr, hits)
    # Separator line
    println(io, SEPARATOR)

    # FBH5 name line
    print(io, "# File ID: ")
    println(io, hdr.fbh5name)
    
    # Separator line
    println(io, SEPARATOR)

    # Source line
    print(io, "# Source: ")
    println(io, hdr.source_name)

    # MJD/RA/Dec line
    print(io, "# MJD: ", hdr.tstart, "\t")
    @printf io "RA: %02dh%02dm%06.3fs\t" ha2hms(hdr.src_raj)...
    @printf io "DEC: %02d:%02d:%06.3f\n" deg2dms(hdr.src_dej)...

    # DELTAF/DELTAT/MAX_DR/OBSLENGTH line
    print(io, "# DELTAT: ", hdr.tsamp, "\t")
    print(io, "DELTAF(Hz): ", hdr.foff*1e6, "\t") # Write as Hz
    print(io, "max_drift_rate: ", hdr.maxdr, "\t")
    println(io, "obs_length: ", hdr.tsamp*hdr.nsamps)

    # Header lines
    println(io, "# --------------------------")
    print(io, "# Top_Hit_#\tDrift_Rate\tSNR\tUncorrected_Frequency\t")
    print(io, "Corrected_Frequency\tIndex\tfreq_start\tfreq_end\t")
    println(io, "SEFD\tSEFD_freq\tCoarse_Channel_Number\tFull_number_of_hits")
    println(io, "# --------------------------")

    # Data table
    for hit in hits
        @printf io "%06d\t" hit.Top_Hit_Number
        @printf io "%.6f\t" hit.Drift_Rate
        @printf io "%.6f\t" hit.SNR
        @printf io "%.6f\t" hit.Uncorrected_Frequency
        @printf io "%.6f\t" hit.Corrected_Frequency
        @printf io "%d\t"   (hit.Index-1) # Write as 0-based
        @printf io "%.6f\t" hit.freq_start
        @printf io "%.6f\t" hit.freq_end
        @printf io "%.1f\t" hit.SEFD
        @printf io "%.6f\t" hit.SEFD_freq
        @printf io "%d\t"   (hit.Coarse_Channel_Number-1) # Write as 0-based
        @printf io "%d"     hit.Full_number_of_hits
        println(io)
    end

    nothing
end

function writedat(datname::AbstractString, hdr, hits)
    # Output DAT file
    open(datname, "w") do io
        writedat(io, hdr, hits)
    end
end

function writedat(dest, hdrhits)
    hdr, hits = hdrhits
    writedat(dest, hdr, hits)
end

@setup_workload begin
    # Create a header and some hits
    hdr = DatFileHeader(
        "dummy.h5", # fbh5name::String
        "Dummy",    # source_name::String
        51544.0,    # tstart::Float64
        0.0,        # src_raj::Float64
        0.0,        # src_dej::Float64
        1.0,        # tsamp::Float64
        1.0,        # foff::Float64
        1.0,        # maxdr::Float64
        10          # nsamps::Int
    )

    # Create some hits
    hits = [
        DatFileHit(1, 0.0, 0.0, 1000.0, 1000.0, 1, 1000.0, 1000.0, 0.0, 0.0, 1, 1)
        DatFileHit(2, 0.0, 0.0, 1000.0, 1000.0, 1, 1000.0, 1000.0, 0.0, 0.0, 1, 1)
    ]

    datname = tempname(; cleanup=false) * ".dat"

    @compile_workload begin
        writedat(datname, (hdr, hits))
        readdat(datname)
        readdat(datname; fix_freq_end=true)
    end
end

end # module SETIDatFiles
