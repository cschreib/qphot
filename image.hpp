struct image_t {
    // Read from parameter file
    std::string filename;
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

    // Work variables
    double lambda = dnan;
    uint_t eazy_band = npos;
    std::string base_filename;
    std::string data_hash;
    std::unique_ptr<fits::input_image> input_image;
    vec2d data;
    vec2d psf;
    fits::header hdr;
    astro::wcs wcs;
    astro::wcs psf_wcs;
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
    double aspix = 0.0;
    double regridded_aspix = 0.0;
    double conv_radius = 0.0;
    double noise_correl = 1.0;
    double flux = dnan, flux_err = dnan, apcor = dnan, background = dnan;
    vec1d flux_bg;
    uint_t num_bg = 0;

    // Helper functions
    void regrid_data(const astro::wcs& det_wcs, const astro::wcs& det_psf_wcs, const fits::header& det_hdr,
        const fits::header& det_psf_hdr, double det_aspix, bool reuse, const std::string& det_hash,
        const std::string& out_dir) {

        // See if we can reuse a cached image
        this->regridded_data_cache = out_dir+this->base_filename+"_regrid_cache.fits";

        std::string cache_hash = hash(this->data_hash, det_hash);

        bool remake_regrid = true;
        if (reuse && file::exists(this->regridded_data_cache)) {
            fits::input_image iimg(this->regridded_data_cache);
            iimg.reach_hdu(1);

            std::string thash;
            if (iimg.read_keyword("CACHEID", thash) && thash == cache_hash) {
                // We're good, read the cached image
                iimg.read(this->regridded_data);
                iimg.reach_hdu(2);
                iimg.read(this->regridded_psf);

                this->regridded_aspix = det_aspix;
                remake_regrid = false;
            }
        }

        if (remake_regrid) {
            // Regrid
            if (det_aspix < aspix) {
                regrid_interpolate_params rp;
                rp.conserve_flux = true;

                // Estimate noise correlation resulting from regriding
                if (rp.method == interpolation_method::linear) {
                    // Empirically estimated with simulations
                    noise_correl = 1.5*sqr(aspix/det_aspix) + 0.78*pow(aspix/det_aspix, -2.3);
                } else if (rp.method == interpolation_method::nearest) {
                    // 1 pixel -> n pixels with the same value, simple
                    noise_correl = sqr(aspix/det_aspix);
                }

                this->regridded_data = regrid_interpolate(this->data, this->wcs, det_wcs, rp);
                // Regrid the PSF as well
                this->regridded_psf = regrid_interpolate(this->psf, this->psf_wcs, det_psf_wcs, rp);
            } else {
                this->regridded_data = regrid_drizzle(this->data, this->wcs, det_wcs);
                // Regrid the PSF as well
                this->regridded_psf = regrid_drizzle(this->psf, this->psf_wcs, det_psf_wcs);
            }

            this->regridded_psf[where(!is_finite(this->regridded_psf))] = 0;
            this->regridded_aspix = det_aspix;

            // Make sure the PSF is still centered
            vec1i idm = mult_ids(this->regridded_psf, max_id(this->regridded_psf));
            this->regridded_psf = recenter(this->regridded_psf, idm[0], idm[1]);

            // Cache them for future use
            fits::output_image oimg(this->regridded_data_cache);
            oimg.reach_hdu(1);
            oimg.write(this->regridded_data);
            oimg.write_header(det_hdr);
            oimg.write_keyword("CACHEID", cache_hash);
            oimg.reach_hdu(2);
            oimg.write(this->regridded_psf);
            oimg.write_header(det_psf_hdr);
        }

        this->idf = where(is_finite(this->regridded_data));
        this->idnf = where(!is_finite(this->regridded_data));

        this->regridded_hdr = det_hdr;

        if (this->hri) {
            this->regridded_data_noconvol = this->regridded_data;
        }

        fits::setkey(regridded_hdr, "NCORR", noise_correl, "[noise correlation]");
    }

    void convolve_to(double new_seeing) {
        if (new_seeing <= seeing) return;

        conv_radius = sqrt(sqr(new_seeing) - sqr(seeing))/2.355/regridded_aspix;

        // Build kernel
        double ky0 = regridded_data.dims[0]/2;
        double kx0 = regridded_data.dims[1]/2;
        vec2d kernel = gaussian_profile(regridded_data.dims, conv_radius, ky0, kx0);
        noise_correl = sqrt(sqrt(noise_correl) + 1.0/total(sqr(kernel)));

        // Convolve image
        make_finite();
        regridded_data = convolve2d(regridded_data, kernel);
        restore_not_finite();

        seeing = new_seeing;

        // Convolve the PSF as well
        regridded_psf = convolve2d(regridded_psf, kernel);

        convolved = true;

        fits::setkey(regridded_hdr, "CONVOL", conv_radius, "[pixels]");
        fits::setkey(regridded_hdr, "NCORR", noise_correl, "[noise correlation]");
    }

    void fit_neighbors(linfit_batch_t<vec1d>& batch, const vec2d& models,
        const vec1u& idfit, const vec1u& gmodel, double dance_chi2_range, int_t dy, int_t dx) {

        vec2d timg;
        if (dy != 0 || dx != 0) {
            // Move image in the opposite direction
            timg = translate_integer(regridded_data, -dy, -dx);
        } else {
            timg = regridded_data;
        }

        // Fit
        batch.fit_nochi2(timg[idfit]);

        // Build local residual and chi2
        // Subtract models
        for (uint_t i : gmodel)
        for (uint_t k : idfit) {
            timg.safe[k] -= batch.fr.params.safe[i]*models.safe(i,k);
        }

        // Comput local chi2
        batch.fr.chi2 = 0; {
            int_t nlpix = ceil(dance_chi2_range/regridded_aspix);
            int_t iyl0 = int_t(timg.dims[0]/2) + dy;
            int_t ixl0 = int_t(timg.dims[1]/2) + dx;
            uint_t yl0 = (iyl0 > nlpix ? iyl0-nlpix : 0);
            uint_t xl0 = (ixl0 > nlpix ? ixl0-nlpix : 0);
            uint_t yl1 = (iyl0 < int_t(timg.dims[0])-nlpix ? iyl0+nlpix : int_t(timg.dims[0]-1));
            uint_t xl1 = (ixl0 < int_t(timg.dims[1])-nlpix ? ixl0+nlpix : int_t(timg.dims[1]-1));

            vec2b masked = replicate(true, timg.dims);
            masked[idfit] = false;
            for (uint_t iy : range(yl1-yl0+1))
            for (uint_t ix : range(xl1-xl0+1)) {
                if (!masked.safe(iy+yl0,ix+xl0)) {
                    batch.fr.chi2 += sqr(timg.safe(iy+yl0,ix+xl0));
                }
            }
        }
    }

    void make_finite() {
        regridded_data[idnf] = 0;
    }

    void restore_not_finite() {
        regridded_data[idnf] = dnan;
    }
};

bool read_image_list(const std::string& filename, vec<1,image_t>& imgs) {
    std::string idir = file::get_directory(filename);

    std::string line;
    std::ifstream file(filename);
    uint_t l = 0;
    image_t* cimg = nullptr;

    while (std::getline(file, line)) {
        ++l;
        line = trim(line);
        if (line.empty() || line[0] == '#') continue;

        auto eqpos = line.find('=');
        if (eqpos == line.npos) {
            std::string filename;
            if (line[0] == '/') {
                filename = line;
            } else {
                filename = idir+line;
            }

            auto iter = std::find_if(imgs.begin(), imgs.end(), [&](const image_t& img) {
                return img.filename == filename;
            });

            if (iter == imgs.end()) {
                image_t img;
                img.filename = filename;
                img.base_filename = file::remove_extension(file::get_basename(img.filename));
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
            } else if (key == "seeing") {
                if (!from_string(val, cimg->seeing)) {
                    write_error("could not read seeing value from '", val, "'");
                    write_error("at ", filename, ":", l);
                    return false;
                }
            } else if (key == "zero_point") {
                if (!from_string(val, cimg->zero_point)) {
                    write_error("could not read zero point value from '", val, "'");
                    write_error("at ", filename, ":", l);
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
                    write_error("could not read 'aperture' value from '", val, "'");
                    write_error("at ", filename, ":", l);
                    return false;
                }
            } else if (key == "nodetect") {
                if (!from_string(val, cimg->nodetect)) {
                    write_error("could not read 'nodetect' value from '", val, "'");
                    write_error("at ", filename, ":", l);
                    return false;
                }
            } else if (key == "clean_neighbors") {
                if (!from_string(val, cimg->clean_neighbors)) {
                    write_error("could not read 'clean_neighbors' value from '", val, "'");
                    write_error("at ", filename, ":", l);
                    return false;
                }
            } else if (key == "noconvolve") {
                if (!from_string(val, cimg->noconvolve)) {
                    write_error("could not read 'noconvolve' value from '", val, "'");
                    write_error("at ", filename, ":", l);
                    return false;
                }
            } else if (key == "nodance") {
                if (!from_string(val, cimg->nodance)) {
                    write_error("could not read 'nodance' value from '", val, "'");
                    write_error("at ", filename, ":", l);
                    return false;
                }
            } else if (key == "flux_method") {
                if (val != "aper" && val != "none") {
                    write_error("unknown flux method '", val, "'");
                    write_error("allowed values are: 'aper', 'none'");
                    return false;
                }
                cimg->flux_method = val;
            } else if (key == "hri") {
                if (!from_string(val, cimg->hri)) {
                    write_error("could not read 'hri' value from '", val, "'");
                    write_error("at ", filename, ":", l);
                    return false;
                }
            } else {
                write_warning("unknown parameter '", key, "', ignored");
                write_warning("at ", filename, ":", l);
            }
        }
    }

    return true;
}
