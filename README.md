# Echo-top Heights Calculation

This project calculates the echo-top heights of radar volumes at
multiple reflectivity thresholds. It is implemented in C++ and utilizes
the BOM (Bureau of Meteorology) library.

## Features

- Read radar data from HDF5 files.
- Compute echo-top heights for various reflectivity thresholds.
- Output results to new radar volume files. 

## Requirements

- C++17 or higher
- BOM libraries
- HDF5 library
- CMake (for building the project)

## Usage

``` sh
echotop [options] vol.h5
```

## Options

- `-h, --help` Show this message and exit.
- `-r, --refl` Reflectivity threshold for the echo-top heights. Can be a
  single value or a list of values separated by a comma (e.g., -r
  1,5.5,10).

## Code Structure

### Headers

The project includes the following headers:

- <span class="title-ref">bom/io/odim.h</span>
- <span class="title-ref">bom/radar/beam_propagation.h</span>
- <span class="title-ref">bom/array2.h</span>
- <span class="title-ref">bom/ellipsoid.h</span>
- <span class="title-ref">bom/grid_coordinates.h</span>
- <span class="title-ref">bom/grid_transform.h</span>
- <span class="title-ref">bom/map_projection.h</span>
- <span class="title-ref">bom/trace.h</span>

### Libraries

The following libraries are used:

- <span class="title-ref">\<algorithm\></span>
- <span class="title-ref">\<filesystem\></span>
- <span class="title-ref">\<getopt.h\></span>
- <span class="title-ref">\<math.h\></span>
- <span class="title-ref">\<numeric\></span>
- <span class="title-ref">\<string\></span>
- <span class="title-ref">\<thread\></span>
- <span class="title-ref">\<vector\></span>

## Functions

### Main Functions

- **get_elevation(io::odim::polar_volume const vol_odim):** Extracts
  elevation angles from the radar volume.
- **read_moment(io::odim::polar_volume const vol_odim, string moment):**
  Reads radar data for a specific moment (e.g., DBZH).
- **read_volume(const std::filesystem::path& path, const string
  reflname):** Reads radar volume data from an HDF5 file.
- **run_search(const array1\<bin_info\>& bins, float target):** Binary
  search for specific bin information.
- **compute_eth(radarset dset, float dbz_thld):** Computes echo-top
  heights for a given reflectivity threshold.
- **process_file(const std::filesystem::path& path, std::vector\<float\>
  r_values):** Processes radar volume data and computes echo-top heights
  for specified reflectivity thresholds.

### Utility Functions

- **argsort(const std::vector\<T\> &array):** Sorts indices based on the
  array values.
- **create_eth_label(float value):** Creates a label for echo-top
  heights based on the reflectivity threshold value.
- **parse_floats(const std::string &s, std::vector\<float\> &v):**
  Parses a string of float values separated by commas.

## Building

To build the project, use CMake:

``` sh
mkdir build
cd build
cmake ..
make
```

## Running the Program

To run the program, use the following command:

``` sh
./echotop [options] vol.h5
```
