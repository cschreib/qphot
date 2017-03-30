void state_t::extract_cutouts(double ra, double dec) {
    if (!good) return;

    if (opts.verbose) {
        note("extracting cutouts...");
    }

    for (auto& img : images) {
        cutout_extractor ex;
        ex.setup_image(img.source.filename);
        if (!img.source.distortfile.empty()) {
            ex.setup_distortion(img.source.distortfile);
        }

        ex.get_cutout(img.data, img.hdr, ra, dec, opts.cutout_size);
        img.wcs = astro::wcs(img.hdr);

        if (!opts.nocuts) {
            // Save cutout on disk
            fits::write(outdir+img.source.short_name+".fits", img.data, img.hdr);
        }
    }
}

void state_t::read_cutouts() {
    if (!good) return;

    if (opts.verbose) {
        note("reading images in memory...");
    }

    for (auto& img : images) {
        fits::input_image iimg(img.source.input_images[0]);
        iimg.read(img.data);
        img.hdr = iimg.read_header();
        img.wcs = astro::wcs(img.hdr);
    }
}

// Measure the "seeing" (FWHM) from a PSF image
double get_seeing(vec2d img, double y0, double x0) {
    auto get_seeing_raw = [](const vec2d& timg, double ty0, double tx0) {
        double tot = 0.0;
        for (uint_t y : range(timg.dims[0]))
        for (uint_t x : range(timg.dims[1])) {
            tot += timg.safe(y,x)*(sqr(y - ty0) + sqr(x - tx0));
        }

        return sqrt(tot);
    };

    // Compute first estimate
    img /= total(img); // image needs to be normalized first
    double rr = get_seeing_raw(img, y0, x0);

    // Simulate various profiles
    vec1d td = rgen_step(0.1, max(img.dims[0], img.dims[1]), 0.3);
    vec1d md(td.dims);

    vec2d timg(img.dims);
    for (uint_t i : range(td)) {
        for (uint_t y : range(timg.dims[0]))
        for (uint_t x : range(timg.dims[1])) {
            timg.safe(y,x) = integrate_gauss_2d(y-0.5, y+0.5, x-0.5, x+0.5, y0, x0, td[i]);
        }

        md[i] = get_seeing_raw(timg, y0, x0);
    }

    // Correct measured value using simulation
    // The additional factor of 2*sqrt(log(4.0)) is the conversion
    // from a Gaussian sigma into a FWHM (=2.355)
    return interpolate(td, md, rr)*2*sqrt(log(4.0));
}

void state_t::read_psfs() {
    if (!good) return;

    for (auto& img : images) {
        auto& idims = img.data.dims;

        if (!img.source.psffile.empty()) {
            // Read the PSF
            fits::read(img.source.psffile, img.psf);

            // Make sure the PSF is centered
            vec1i idm = mult_ids(img.psf, max_id(img.psf));
            img.psf = recenter(img.psf, idm[0], idm[1]);

            if (!is_finite(img.seeing)) {
                img.seeing = get_seeing(img.psf,
                    img.psf.dims[0]/2, img.psf.dims[1]/2);
            }
        } else {
            // Create a PSF if it is not provided
            // NB: make sure it is large enough that any source in the image will be
            // properly modeled without reaching the edge of the PSF image
            img.psf = gaussian_profile({{idims[0]*2+1, idims[1]*2+1}},
                img.seeing/2.355/img.source.aspix, idims[0], idims[1]
            );
        }

        {
            // Create astrometry of the PSF image
            double tra, tdec;
            astro::xy2ad(img.wcs, idims[1]/2+1, idims[0]/2+1, tra, tdec);

            fits::header thdr = img.hdr;
            fits::setkey(thdr, "NAXIS1", img.psf.dims[1]);
            fits::setkey(thdr, "NAXIS2", img.psf.dims[0]);
            fits::setkey(thdr, "CRPIX1", img.psf.dims[1]/2 + 1);
            fits::setkey(thdr, "CRPIX2", img.psf.dims[0]/2 + 1);
            fits::setkey(thdr, "CRVAL1", tra);
            fits::setkey(thdr, "CRVAL2", tdec);
            img.psf_wcs = astro::wcs(thdr);
        }
    }
}
