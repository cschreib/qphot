void state_t::regrid(image_t& img) {
    // See if we can reuse a cached image
    img.regridded_data_cache = outdir+img.source.short_name+"_regrid_cache.fits";

    std::string cache_hash = hash(img.data_hash, det_cache_hash);

    bool remake_regrid = true;
    if (opts.reuse && file::exists(img.regridded_data_cache)) {
        fits::input_image iimg(img.regridded_data_cache);
        iimg.reach_hdu(1);

        std::string thash;
        if (iimg.read_keyword("CACHEID", thash) && thash == cache_hash) {
            // We're good, read the cached image
            iimg.read(img.regridded_data);
            iimg.reach_hdu(2);
            iimg.read(img.regridded_psf);

            img.regridded_aspix = det_aspix;
            remake_regrid = false;
        }
    }

    if (remake_regrid) {
        // Regrid
        if (det_aspix < img.source.aspix) {
            regrid_interpolate_params rp;
            rp.conserve_flux = true;

            // Estimate noise correlation resulting from regriding
            if (rp.method == interpolation_method::linear) {
                // Empirically estimated with simulations
                img.noise_correl = 1.5*sqr(img.source.aspix/det_aspix) + 0.78*pow(img.source.aspix/det_aspix, -2.3);
            } else if (rp.method == interpolation_method::nearest) {
                // 1 pixel -> n pixels with the same value, simple
                img.noise_correl = sqr(img.source.aspix/det_aspix);
            }

            img.regridded_data = regrid_interpolate(img.data, img.wcs, det_wcs, rp);
            // Regrid the PSF as well
            img.regridded_psf = regrid_interpolate(img.psf, img.psf_wcs, det_psf_wcs, rp);
        } else {
            img.regridded_data = regrid_drizzle(img.data, img.wcs, det_wcs);
            // Regrid the PSF as well
            img.regridded_psf = regrid_drizzle(img.psf, img.psf_wcs, det_psf_wcs);
        }

        img.regridded_psf[where(!is_finite(img.regridded_psf))] = 0;
        img.regridded_aspix = det_aspix;

        // Make sure the PSF is still centered
        vec1i idm = mult_ids(img.regridded_psf, max_id(img.regridded_psf));
        img.regridded_psf = recenter(img.regridded_psf, idm[0], idm[1]);

        // Cache them for future use
        if (!opts.nocuts) {
            fits::output_image oimg(img.regridded_data_cache);
            oimg.reach_hdu(1);
            oimg.write(img.regridded_data);
            oimg.write_header(det_hdr);
            oimg.write_keyword("CACHEID", cache_hash);
            oimg.reach_hdu(2);
            oimg.write(img.regridded_psf);
            oimg.write_header(det_psf_hdr);
        }
    }

    img.idf = where(is_finite(img.regridded_data));
    img.idnf = where(!is_finite(img.regridded_data));

    img.regridded_hdr = det_hdr;

    if (img.source.hri && (opts.clean_model == "hri" || opts.clean_model == "nearest_hri")) {
        img.regridded_data_noconvol = img.regridded_data;
    }

    fits::setkey(img.regridded_hdr, "NCORR", img.noise_correl, "[noise correlation]");
}

void state_t::convolve_to(image_t& img, double new_seeing) {
    if (new_seeing <= img.seeing) {
        img.convolved = true;
        return;
    }

    img.conv_radius = sqrt(sqr(new_seeing) - sqr(img.seeing))/2.355/img.regridded_aspix;

    // Build kernel
    double ky0 = img.regridded_data.dims[0]/2;
    double kx0 = img.regridded_data.dims[1]/2;
    vec2d kernel = gaussian_profile(img.regridded_data.dims, img.conv_radius, ky0, kx0);
    img.noise_correl = sqrt(sqrt(img.noise_correl) + 1.0/total(sqr(kernel)));

    // Convolve image
    make_finite(img);
    img.regridded_data = convolve2d(img.regridded_data, kernel);
    restore_not_finite(img);
    img.seeing = new_seeing;

    // Convolve the PSF as well
    img.regridded_psf = convolve2d(img.regridded_psf, kernel);

    img.convolved = true;

    fits::setkey(img.regridded_hdr, "CONVOL", img.conv_radius, "[pixels]");
    fits::setkey(img.regridded_hdr, "NCORR", img.noise_correl, "[noise correlation]");
}

void state_t::fit_neighbors(image_t& img, linfit_batch_t<vec1d>& batch, const vec2d& models,
    const vec1u& idfit, const vec1u& gmodel, int_t dy, int_t dx) const {

    vec2d timg;
    if (dy != 0 || dx != 0) {
        // Move image in the opposite direction
        timg = translate_integer(img.regridded_data, -dy, -dx);
    } else {
        timg = img.regridded_data;
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
        int_t nlpix = ceil(opts.dance_chi2_range/img.regridded_aspix);
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

void state_t::make_finite(image_t& img) {
    img.regridded_data[img.idnf] = 0;
}

void state_t::restore_not_finite(image_t& img) {
    img.regridded_data[img.idnf] = dnan;
}
