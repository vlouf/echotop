#include <bom/io/configuration.h>
#include <bom/io/cf.h>
#include <bom/io/nc.h>
#include <bom/io/odim.h>
#include <bom/radar/beam_propagation.h>
#include <bom/array2.h>
#include <bom/ellipsoid.h>
#include <bom/grid_coordinates.h>
#include <bom/grid_transform.h>
#include <bom/map_projection.h>
#include <bom/trace.h>

#include <algorithm>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <getopt.h>
#include <iomanip>
#include <iostream>
#include <math.h>
#include <numeric>
#include <sstream>
#include <string>
#include <thread>
#include <unordered_map>
#include <vector>

using namespace bom;

constexpr float nodata = std::numeric_limits<float>::quiet_NaN();
constexpr float undetect = -32.0f;

struct bin_info{
  float slant_range;
  float ground_range;
  float altitude;
};

struct sweep{
  radar::beam_propagation beam;
  array1<bin_info>        bins; // @ bin centers
  array1<angle>           rays; // @ ray centers
  array2f                 data;
};

struct volume{
  latlonalt     location;
  vector<sweep> sweeps;
};

struct radarset{
    volume dbzh;
    vector<float> elevation;
    string source;
    string date;
    string time;
  };


constexpr auto try_again = "try --help for usage instructions\n";
constexpr auto usage_string =
R"(Echo-top heights of radar volume at multiple reflectivity

usage:
  echotop [options] vol.h5 out.nc

available options:
  -h, --help
      Show this message and exit

  -g, --generate
      Output a sample configuration file and exit

  -t, --trace=level
      Set logging level [log]
        none | status | error | warning | log | debug
)";

constexpr auto short_options = "hgt:";
constexpr struct option long_options[] =
{
    { "help",     no_argument,       0, 'h' }
  , { "generate", no_argument,       0, 'g' }
  , { "trace",    required_argument, 0, 't' }
  , { 0, 0, 0, 0 }
};

template<typename T>
std::vector<size_t> argsort(const std::vector<T> &array) {
    std::vector<size_t> indices(array.size());
    std::iota(indices.begin(), indices.end(), 0);
    std::sort(indices.begin(), indices.end(),
              [&array](int left, int right) -> bool {
                  // sort indices according to corresponding array element
                  return array[left] < array[right];
              });

    return indices;
}

auto get_elevation(io::odim::polar_volume const vol_odim) -> vector<float>{
    const size_t nelev = vol_odim.scan_count();
    vector<float> elevation;

    for (size_t iscan = 0; iscan < nelev; ++iscan){
          auto scan_odim = vol_odim.scan_open(iscan);
          elevation.push_back(scan_odim.elevation_angle());
      }

    return elevation;
  }

auto read_moment(io::odim::polar_volume const vol_odim, string moment) -> volume{
    auto vol = volume{};

    vol.location.lat = vol_odim.latitude() * 1_deg;
    vol.location.lon = vol_odim.longitude() * 1_deg;
    vol.location.alt = vol_odim.height();

    for (size_t iscan = 0; iscan < vol_odim.scan_count(); ++iscan)
    {
      auto scan_odim = vol_odim.scan_open(iscan);
      auto scan = sweep{};

      scan.beam = radar::beam_propagation{vol.location.alt, scan_odim.elevation_angle() * 1_deg};
      scan.bins.resize(scan_odim.bin_count());

      auto range_scale = scan_odim.range_scale();
      auto range_start = scan_odim.range_start() * 1000 + range_scale * 0.5;
      for (size_t i = 0; i < scan.bins.size(); ++i)
      {
        scan.bins[i].slant_range = range_start + i * range_scale;
        std::tie(scan.bins[i].ground_range, scan.bins[i].altitude) = scan.beam.ground_range_altitude(scan.bins[i].slant_range);
      }

      scan.rays.resize(scan_odim.ray_count());
      auto ray_scale = 360_deg / scan.rays.size();
      auto ray_start = scan_odim.ray_start() * 1_deg + ray_scale * 0.5;
      for (size_t i = 0; i < scan.rays.size(); ++i)
        scan.rays[i] = ray_start + i * ray_scale;

      for (size_t idata = 0; idata < scan_odim.data_count(); ++idata)
      {
        auto data_odim = scan_odim.data_open(idata);
        if (data_odim.quantity() != moment)
          continue;

        scan.data.resize(vec2z{(size_t) scan_odim.bin_count(), (size_t) scan_odim.ray_count()});
        data_odim.read_unpack(scan.data.data(), undetect, nodata);
        vol.sweeps.push_back(std::move(scan));
        break;
      }
    }

    return vol;
  }

auto read_volume(const std::filesystem::path& path, const string reflname) -> radarset {
    radarset dset;
    io::odim::polar_volume vol_odim{path, io_mode::read_only};

    std::cout << "Reading file " << path << "\n";
    dset.elevation = get_elevation(vol_odim);
    std::cout << "Reading Moment " << reflname << "\n";
    dset.dbzh = read_moment(vol_odim, reflname);

    const auto& attributes = vol_odim.attributes();
    dset.source = attributes["source"].get_string();
    dset.date = attributes["date"].get_string();
    dset.time = attributes["time"].get_string();

    return dset;
}

auto run_search(const array1<bin_info>& bins, float target) -> size_t{
    size_t left = 0;
    size_t right = bins.size() - 2;
    std::cout << "Bin size: " << right << std::endl;

    while (left < right) {
        size_t mid = (right + left) / 2;
        if (bins[mid + 1].ground_range < target) {
            left = mid + 1;
        } else {
            right = mid;
        }
    }
    return left;
}


void process_file(const std::filesystem::path& path){
    float eth_thld = 10.;
    float noise_thld = -5.f;

    std::cout << "Starting process." << std::endl;
    auto dset = read_volume(path, "DBZH");
    auto sorted_index = argsort(dset.elevation);
    const auto nbins = dset.dbzh.sweeps[sorted_index[0]].bins.size();
    const auto nrays = dset.dbzh.sweeps[sorted_index[0]].rays.size();
    auto eth = array2f(vec2z(nbins, nrays));
    eth.fill(0.f);

    auto range_base = dset.dbzh.sweeps[sorted_index[0]].bins;

    for(auto k = 0; k < sorted_index.size() - 1; k ++){
        auto thetab = dset.elevation[sorted_index[k]];  // Ref - following Laks conventions.
        auto thetaa = dset.elevation[sorted_index[k+1]];  // Iter
        auto range_ref = dset.dbzh.sweeps[sorted_index[k]].bins;
        auto range_iter = dset.dbzh.sweeps[sorted_index[k + 1]].bins;
        auto refl_ref = dset.dbzh.sweeps[sorted_index[k]].data;
        auto refl_iter = dset.dbzh.sweeps[sorted_index[k + 1]].data;

        for(auto i = 0; i < nbins; i ++){            
            auto i_r = run_search(range_ref, range_base[i].ground_range);
            auto i_i = run_search(range_iter, range_base[i].ground_range);
            for(auto j = 0; j < nrays; j ++){
                auto refb = refl_ref[j][i_r];
                auto refa = refl_iter[j][i_i];
                if(std::isnan(refa)) refa = noise_thld;
                if(refb >= eth_thld){
                    eth[j][i] = thetab + 0.5;
                    continue;
                }
                if((refb > eth_thld) && (refa <= eth_thld)){
                    if(std::fabs(90. - thetaa ) > 0.1){
                        eth[j][i] = (eth_thld - refa) * (thetab - thetaa) / (refb - refa) + thetab;
                    } else {
                        eth[j][i] = thetab + 0.5;
                    }
                }
            }
        }
    }

    size_t dims[2] = {nrays, nbins};
    io::odim::polar_volume vol_odim{path, io_mode::read_write};
    auto scan_odim = vol_odim.scan_open(sorted_index[0]);
    auto data = scan_odim.data_append(
        io::odim::data::data_type::f32
      , 2
      , dims      
      );
    data.write(eth.data());
    data.set_quantity("eth");
    data.set_nodata(-1.);
    data.set_undetect(-1.);
    data.set_gain(1);
    data.set_offset(0);

    std::cout << "Process done." << std::endl;
}

int main(int argc, char* argv[])
{
  try
  {
    // process command line
    while (true)
    {
      int option_index = 0;
      int c = getopt_long(argc, argv, short_options, long_options, &option_index);
      if (c == -1)
        break;
      switch (c)
      {
      case 'h':
        std::cout << usage_string;
        return EXIT_SUCCESS;
      case 'g':
        std::cout << usage_string;
        return EXIT_SUCCESS;
      case 't':
        trace::set_min_level(from_string<trace::level>(optarg));
        break;
      case '?':
        std::cerr << try_again;
        return EXIT_FAILURE;
      }
    }

    if (argc - optind != 1)
    {
      std::cerr << "missing required parameter\n" << try_again;
      return EXIT_FAILURE;
    }

    process_file(argv[optind]);
  }
  catch (std::exception& err)
  {
    trace::error("fatal exception: {}", format_exception(err));
    return EXIT_FAILURE;
  }
  catch (...)
  {
    trace::error("fatal exception: (unknown exception)");
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}



