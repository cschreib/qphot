#include "image.hpp"
#include "ds9.hpp"

struct options_t {
    std::string outdir = "./";
    std::string extra;
    bool verbose = false;
    bool show_progress = false;
    double aperture = dnan; // [arcsec]
    double search_radius = dnan; // [arcsec]
    std::string exclude;
    std::string include;
    double det_seeing = 0;
    uint_t det_min_area = 30;
    double det_threshold = 0.5;
    double bg_threshold = 2.0;
    std::string error_file;
    uint_t min_bg_aper = 10;
    bool no_neighbor_mask_background = false;
    bool reuse = false;
    bool reuse_cutouts = false;
    std::string clean_model = "clean";
    double clean_ratio = 0.1;
    double clean_threshold = 3.0;
    double clean_conv_seeing = dnan;
    double hri_min_snr = 5.0;
    bool save_models = false;
    bool debug_clean = false;
    bool clean_group = false;
    bool model_dance = false;
    double dance_step = 0.1;
    double dance_max = 0.5;
    double dance_chi2_range = 6.0;
    bool nocuts = false;
    std::string priors;
    double cutout_size = 20.0;
    uint_t save_interval = 0;
    uint_t threads = 1;
    bool merge = false;
    bool update = false;
    double prior_ra = dnan;
    double prior_dec = dnan;
    vec1s only_sources;
};

struct image_t {
    image_t(const image_source_t& s) : source(s) {
        seeing = source.seeing;
        aperture = source.aperture;
    }

    // Source data
    const image_source_t& source;

    // Work variables
    fits::header hdr;
    astro::wcs wcs;
    astro::wcs psf_wcs;
    std::string data_hash;
    vec2d data;
    vec2d psf;
    std::string regridded_data_cache;
    vec2d regridded_data;
    vec2d regridded_data_noconvol;
    vec1u idf;
    vec1u idnf;
    vec2d regridded_psf;
    fits::header regridded_hdr;
    vec2d mask;
    bool covered = true;
    bool convolved = false;
    bool cleaned = false;
    double seeing = dnan;
    double aperture = dnan;
    double regridded_aspix = 0.0;
    double conv_radius = 0.0;
    double noise_correl = 1.0;
    double flux = dnan, flux_err = dnan, apcor = dnan, background = dnan;
    vec1d flux_bg;
    uint_t num_bg = 0;
};

template<typename T>
void swallow_arguments(std::initializer_list<T>) {}

struct state_t {
    const options_t& opts;
    vec<1,image_t> images;
    bool good = true;
    mutable std::ofstream errfile;

    // Local config (overrides 'opts')
    std::string outdir;
    double det_seeing = dnan;
    double aperture = dnan;
    double search_radius = dnan;

    // Detection stage output
    double smallest_aspix = dinf;
    double worst_seeing = 0;
    std::vector<image_t*> hri;
    std::string det_cache_hash;
    vec2d det_img;
    vec2d det_psf;
    astro::wcs det_wcs;
    astro::wcs det_psf_wcs;
    double det_aspix = 0;
    fits::header det_hdr;
    fits::header det_psf_hdr;
    bool has_clean = false;
    bool has_homogenize = false;
    bool has_filters = false;
    uint_t isource = npos;
    vec2d det_model;
    segment_deblend_output segments;
    vec2u seg;

    state_t(const vec<1,image_source_t>& sources, const options_t& o) : opts(o) {
        images.data.reserve(sources.size());
        for (auto& simg : sources) {
            images.data.emplace_back(simg);
        }

        images.dims[0] = sources.size();

        det_seeing = o.det_seeing;
        aperture = o.aperture;
        search_radius = o.search_radius;
        outdir = o.outdir;
    }

    template<typename T>
    void write_impl(const T& msg) const {
        errfile << msg;
    }

    template<typename ... Args>
    void write_errfile(const Args& ... args) const {
        if (!errfile.is_open()) return;

        // Print time stamp
        std::time_t t = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
        std::tm tm = *std::localtime(&t);
        errfile << "["
                << align_right(strn(tm.tm_year+1900),4,'0') << "/"
                << align_right(strn(tm.tm_mon+1),2,'0') << "/"
                << align_right(strn(tm.tm_mday),2,'0') << " "
                << align_right(strn(tm.tm_hour),2,'0') << ":"
                << align_right(strn(tm.tm_min),2,'0') << ":"
                << align_right(strn(tm.tm_sec),2,'0') << "]: ";

        // Print arguments
        swallow_arguments({(write_impl(args), 0)...});

        // End line
        errfile << std::endl;
    }

    template<typename ... Args>
    void write_error(const Args& ... args) const {
        if (opts.verbose) error(args...);
        write_errfile("error: ", args...);
    }

    template<typename ... Args>
    void write_warning(const Args& ... args) const {
        if (opts.verbose) warning(args...);
        write_errfile("warning: ", args...);
    }

    // state_image.hpp
    void regrid(image_t& img);
    void convolve_to(image_t& img, double new_seeing);
    void fit_neighbors(image_t& img, linfit_batch_t<vec1d>& batch, const vec2d& models,
        const vec1u& idfit, const vec1u& gmodel, int_t dy, int_t dx) const;
    void make_finite(image_t& img);
    void restore_not_finite(image_t& img);
    // state_cutouts.hpp
    void extract_cutouts(double ra, double dec);
    void read_cutouts();
    void read_psfs();
    // state_detect.hpp
    void build_detection_image();
    // state_filters.hpp
    void apply_filters();
    // state_clean.hpp
    void clean_images();
    // state_extract.hpp
    void extract_fluxes();
};

vec2d make_point_source_integer(const vec2d& psf, const std::array<uint_t,2>& dims, int_t y, int_t x) {
    vec2d model(dims);
    int_t y0 = max(0, y - int_t(psf.dims[0]/2));
    int_t y1 = min(int_t(dims[0]), y + int_t(psf.dims[0]/2+1));
    int_t x0 = max(0, x - int_t(psf.dims[1]/2));
    int_t x1 = min(int_t(dims[0]), x + int_t(psf.dims[1]/2+1));

    for (int_t iy = y0; iy < y1; ++iy)
    for (int_t ix = x0; ix < x1; ++ix) {
        int_t ipy = iy - y + int_t(psf.dims[0]/2);
        int_t ipx = ix - x + int_t(psf.dims[1]/2);
        model.safe(uint_t(iy),uint_t(ix)) = psf.safe(uint_t(ipy),uint_t(ipx));
    }

    return model;
}

vec2d make_point_source(const vec2d& psf, const std::array<uint_t,2>& dims, double y, double x) {
    int_t iy = round(y), ix = round(x);
    vec2d npsf = translate(psf, y - iy, x - ix);
    return make_point_source_integer(npsf, dims, iy, ix);
}

#include "state_image.hpp"
#include "state_cutouts.hpp"
#include "state_detect.hpp"
#include "state_filters.hpp"
#include "state_clean.hpp"
#include "state_extract.hpp"
