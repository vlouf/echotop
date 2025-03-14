#include <bom/io/odim.h>
#include <bom/radar/beam_propagation.h>
#include <bom/array2.h>
#include <bom/ellipsoid.h>
#include <bom/grid_coordinates.h>
#include <bom/grid_transform.h>
#include <bom/map_projection.h>
#include <bom/trace.h>

#include <algorithm>
#include <filesystem>
#include <getopt.h>
#include <math.h>
#include <numeric>
#include <string>
#include <thread>
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
  echotop [options] vol.h5

available options:
  -h, --help
      Show this message and exit

  -r, --refl
      Reflectivity threshold for the echo top heights. Can be a single value
      or a list of values separated by a coma (e.g.: -r 1,5.5,10)

  -t, --trace=level
      Set logging level [log]
        none | status | error | warning | log | debug
)";

constexpr auto short_options = "hgt:";
constexpr struct option long_options[] =
{
    { "help",     no_argument,       0, 'h' }
  , { "generate", no_argument,       0, 'g' }
  , { "refl",     required_argument, 0, 'r' }
  , { "trace",    required_argument, 0, 't' }
  , { 0, 0, 0, 0 }
};

template<typename T>
std::vector<size_t> argsort(const std::vector<T> &array) {  // From stackoverflow.
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

  for (size_t iscan = 0; iscan < vol_odim.scan_count(); ++iscan){
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

auto compute_eth(radarset dset, float dbz_thld) -> array2f{
  const double ea = 6371000; // Earth radius in meters.
  float noise_thld = dbz_thld < -2.f ? dbz_thld - 2. : -2.f;
  auto sorted_index = argsort(dset.elevation);
  const auto nbins = dset.dbzh.sweeps[sorted_index[0]].bins.size();
  const auto nrays = dset.dbzh.sweeps[sorted_index[0]].rays.size();
  auto eth = array2f(vec2z(nbins, nrays));
  eth.fill(0.f);

  auto range_base = dset.dbzh.sweeps[sorted_index[0]].bins;

  // Get the apparent elevation angle at which the ETH is found.
  for(size_t k = 0; k < sorted_index.size() - 1; k ++){
    auto thetab = dset.elevation[sorted_index[k]];  // Ref - following Laks conventions.
    auto thetaa = dset.elevation[sorted_index[k+1]];  // Iter
    auto range_ref = dset.dbzh.sweeps[sorted_index[k]].bins;
    auto range_iter = dset.dbzh.sweeps[sorted_index[k + 1]].bins;
    auto refl_ref = dset.dbzh.sweeps[sorted_index[k]].data;
    auto refl_iter = dset.dbzh.sweeps[sorted_index[k + 1]].data;

    for(size_t i = 0; i < range_iter.size(); i ++){
      auto i_b = run_search(range_ref, range_base[i].ground_range);
      auto i_a = run_search(range_iter, range_base[i].ground_range);

      for(size_t j = 0; j < nrays; j ++){
        auto refb = refl_ref[j][i_b];
        auto refa = refl_iter[j][i_a];
        if(std::isnan(refa) || (refa < noise_thld)) refa = noise_thld;

        if((refb > dbz_thld) && (refa < dbz_thld)){
          if(thetaa < 90){
            eth[j][i] = (dbz_thld - refa) * (thetab - thetaa) / (refb - refa) + thetab;
          } else {
            eth[j][i] = thetab + 0.5;
          }
        }
      }
    }
  }

  // Convert the apparent elevation angle into meters MSL and correct for 4/3 Earth radius model.
  for(size_t i = 0; i < nbins; i ++){
    for(size_t j = 0; j < nrays; j ++){
      if(eth[j][i] ==  0) continue;
      double r = range_base[i].ground_range;
      eth[j][i] = sqrt(r * r + 16./9. * ea * ea + 2 * r * 4./3. * ea * sin(M_PI * eth[j][i] / 180.)) - 4./3. * ea;
    }
  }

  return eth;
}

std::string create_eth_label(float value) {
  std::ostringstream oss;
  oss << "ETH_" << value;
  return oss.str();
}

void parse_floats(const std::string &s, std::vector<float> &v) {
  std::stringstream ss(s);
  std::string item;
  while (std::getline(ss, item, ',')) {
      v.push_back(std::stof(item));
  }
}

void process_file(const std::filesystem::path& path, std::vector<float> r_values){
  std::cout << "Starting process." << std::endl;
  auto dset = read_volume(path, "DBZH");
  auto sorted_index = argsort(dset.elevation);
  const auto nbins = dset.dbzh.sweeps[sorted_index[0]].bins.size();
  const auto nrays = dset.dbzh.sweeps[sorted_index[0]].rays.size();
  size_t dims[2] = {nrays, nbins};
  vector<array2f> eth_lst;

  for(float dbz_thld : r_values){
    std::cout << "Computing " << dbz_thld << "dB echo top heights.\n";
    auto eth = compute_eth(dset, dbz_thld);
    eth_lst.push_back(eth);
  }

  io::odim::polar_volume vol_odim{path, io_mode::read_write};
  auto scan_odim = vol_odim.scan_open(sorted_index[0]);
  for(size_t i=0; i<r_values.size(); i++){
    auto dbz_thld = r_values[i];
    auto data = scan_odim.data_append(io::odim::data::data_type::f32, 2, dims);
    data.write(eth_lst[i].data());
    data.set_quantity(create_eth_label(dbz_thld));
    data.set_nodata(-1.);
    data.set_undetect(-1.);
    data.set_gain(1);
    data.set_offset(0);
    std::cout << dbz_thld << "dB ETH written.\n";
  }
  std::cout << "Process done." << std::endl;
}

int main(int argc, char* argv[]){
  try{
    std::vector<float> r_values;
    // process command line
    while (true){
      int option_index = 0;
      int c = getopt_long(argc, argv, short_options, long_options, &option_index);
      if (c == -1)
        break;
      switch (c){
        case 'h':
          std::cout << usage_string;
          return EXIT_SUCCESS;        
        case 'r':
          parse_floats(optarg, r_values);
          break;
        case 't':
          trace::set_min_level(from_string<trace::level>(optarg));
          break;
        case '?':
          std::cerr << try_again;
          return EXIT_FAILURE;
      }
    }

  if (argc - optind != 1){
    std::cerr << "missing required parameter\n" << try_again;
    return EXIT_FAILURE;
  }
  if(r_values.size() == 0) r_values.push_back(0.);  // Default value of 0 dB if none provided.

  process_file(argv[optind], r_values);
  }catch (std::exception& err){
    trace::error("fatal exception: {}", format_exception(err));
    return EXIT_FAILURE;
  }catch (...){
    trace::error("fatal exception: (unknown exception)");
    return EXIT_FAILURE;
  }
}
