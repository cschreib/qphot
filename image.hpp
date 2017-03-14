struct input_image_t {
    // FITS data
    std::unique_ptr<fits::input_image> fits;
    // Cached data
    vec1u dims;
    fits::header hdr;
    astro::wcs wcs;
};

struct image_source_t {
    // Read from parameter file
    std::string filename;
    std::string short_name;
    std::string band;
    std::string psffile;
    double seeing = dnan;
    double zero_point = 23.9;
    bool nodetect = false;
    double aperture = dnan;
    bool clean_neighbors = false;
    bool noconvolve = false;
    bool nodance = false;
    bool hri = false;
    std::string flux_method = "none";

    // Static variables
    double lambda = dnan;
    uint_t eazy_band = npos;
    vec<1,input_image_t> input_images;
    double aspix = dnan;
};

bool read_image_source_list(const std::string& filename, vec<1,image_source_t>& imgs) {
    std::string idir = file::get_directory(filename);

    std::string line;
    std::ifstream file(filename);
    uint_t l = 0;
    image_source_t* cimg = nullptr;

    while (std::getline(file, line)) {
        ++l;
        line = trim(line);
        if (line.empty() || line[0] == '#') continue;

        auto eqpos = line.find('=');
        if (eqpos == line.npos) {
            std::string fname;
            if (line[0] == '/') {
                fname = line;
            } else {
                fname = idir+line;
            }

            auto iter = std::find_if(imgs.begin(), imgs.end(), [&](const image_source_t& img) {
                return img.filename == fname;
            });

            if (iter == imgs.end()) {
                image_source_t img;
                img.filename = fname;
                img.short_name = file::remove_extension(file::get_basename(img.filename));
                imgs.push_back(std::move(img));
                cimg = &imgs.back();
            } else {
                cimg = &*iter;
            }
        } else {
            std::string key = trim(line.substr(0, eqpos));
            std::string val = trim(line.substr(eqpos+1));
            if (key == "band") {
                cimg->band = val;
            } else if (key == "short_name") {
                cimg->short_name = val;
            } else if (key == "seeing") {
                if (!from_string(val, cimg->seeing)) {
                    error("could not read seeing value from '", val, "'");
                    error("at ", filename, ":", l);
                    return false;
                }
            } else if (key == "zero_point") {
                if (!from_string(val, cimg->zero_point)) {
                    error("could not read zero point value from '", val, "'");
                    error("at ", filename, ":", l);
                    return false;
                }
            } else if (key == "psf") {
                if (!val.empty()) {
                    if (val[0] == '/') {
                        cimg->psffile = val;
                    } else {
                        cimg->psffile = idir+val;
                    }
                }
            } else if (key == "aperture") {
                if (!from_string(val, cimg->aperture)) {
                    error("could not read 'aperture' value from '", val, "'");
                    error("at ", filename, ":", l);
                    return false;
                }
            } else if (key == "nodetect") {
                if (!from_string(val, cimg->nodetect)) {
                    error("could not read 'nodetect' value from '", val, "'");
                    error("at ", filename, ":", l);
                    return false;
                }
            } else if (key == "clean_neighbors") {
                if (!from_string(val, cimg->clean_neighbors)) {
                    error("could not read 'clean_neighbors' value from '", val, "'");
                    error("at ", filename, ":", l);
                    return false;
                }
            } else if (key == "noconvolve") {
                if (!from_string(val, cimg->noconvolve)) {
                    error("could not read 'noconvolve' value from '", val, "'");
                    error("at ", filename, ":", l);
                    return false;
                }
            } else if (key == "nodance") {
                if (!from_string(val, cimg->nodance)) {
                    error("could not read 'nodance' value from '", val, "'");
                    error("at ", filename, ":", l);
                    return false;
                }
            } else if (key == "flux_method") {
                if (val != "aper" && val != "none") {
                    error("unknown flux method '", val, "'");
                    error("allowed values are: 'aper', 'none'");
                    return false;
                }
                cimg->flux_method = val;
            } else if (key == "hri") {
                if (!from_string(val, cimg->hri)) {
                    error("could not read 'hri' value from '", val, "'");
                    error("at ", filename, ":", l);
                    return false;
                }
            } else {
                warning("unknown parameter '", key, "', ignored");
                warning("at ", filename, ":", l);
            }
        }
    }

    return true;
}

bool check_image_source_list(vec<1,image_source_t>& imgs) {
    // Read filter data base
    auto fdb = read_filter_db(data_dir+"fits/filter-db/db.dat");
    auto fmap = read_filter_map(data_dir+"fits/filter-db/fast.map");

    // Get basic properties and checks about images
    for (auto& img : imgs) {
        vec1s sects;
        if (end_with(img.filename, ".sectfits")) {
            sects = fits::read_sectfits(img.filename);
        } else {
            sects.push_back(img.filename);
        }

        img.input_images.resize(sects.size());
        for (uint_t i : range(sects)) {
            auto& iimg = img.input_images[i];
            iimg.fits = std::unique_ptr<fits::input_image>(
                new fits::input_image(sects[i])
            );

            if (!iimg.fits->is_image()) {
                error("'", sects[i], "' is not a FITS image");
                if (sects.size() > 1) {
                    error("reading '", img.filename, "'");
                }
                return false;
            }

            iimg.dims = iimg.fits->image_dims();
            iimg.hdr = iimg.fits->read_header();
            iimg.wcs = astro::wcs(iimg.hdr);

            double aspix = dnan;
            if (!astro::get_pixel_size(iimg.wcs, aspix)) {
                error("could not find pixel size in '", sects[i], "'");
                if (sects.size() > 1) {
                    error("reading '", img.filename, "'");
                }
                return false;
            }

            if (!is_finite(img.aspix)) {
                img.aspix = aspix;
            } else if (abs(img.aspix - aspix) > 0.001) {
                error("pixel size is not homogeneous in '", img.filename, "'");
                return false;
            }
        }

        filter_t fil;
        if (!get_filter(fdb, img.band, fil)) {
            warning("unknown filter '", img.band, "'");
        } else {
            uint_t eid = npos;
            get_filter_id(fmap, img.band, eid);

            img.lambda = fil.rlam;
            img.eazy_band = eid;
        }
    }

    return true;
}
