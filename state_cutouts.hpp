void state_t::extract_cutouts(double ra, double dec) {
    if (!good) return;

    if (opts.verbose) {
        note("extracting cutouts...");
    }

    for (auto& img : images) {
        int_t hs = max(10, ceil(0.5*opts.cutout_size/img.source.aspix));
        double px = hs+1.0, py = hs+1.0;

        fits::header first_hdr;
        for (uint_t i : range(img.source.input_images)) {
            // Look on each FITS image of the set if the source is covered
            fits::input_image iimg(img.source.input_images[i]);
            fits::header hdr = iimg.read_header();
            if (first_hdr.empty()) {
                first_hdr = hdr;
            }

            astro::wcs w = astro::wcs(hdr);
            vec1u dims = iimg.image_dims();

            double dxc, dyc;
            astro::ad2xy(w, ra, dec, dxc, dyc);
            dxc -= 1.0; dyc -= 1.0;
            int_t xc = round(dxc);
            int_t yc = round(dyc);
            int_t iy0 = yc-hs, iy1 = yc+hs, ix0 = xc-hs, ix1 = xc+hs;

            if (ix1 < 0 || ix0 >= int_t(dims[1]) || iy1 < 0 || iy0 >= int_t(dims[0])) {
                // Source not covered
                continue;
            }

            vec2d data;
            if (ix0 < 0 || ix1 >= int_t(dims[1]) || iy0 < 0 || iy1 >= int_t(dims[0])) {
                // Source partially covered
                data = replicate(dnan, 2*hs+1, 2*hs+1);

                uint_t tx0 = max(0, ix0);
                uint_t ty0 = max(0, iy0);
                uint_t tx1 = min(dims[1]-1, ix1);
                uint_t ty1 = min(dims[0]-1, iy1);

                vec2d subcut;
                iimg.read_subset(subcut, ty0-_-ty1, tx0-_-tx1);

                uint_t x0 = int_t(tx0)-ix0, x1 = int_t(tx1)-ix0, y0 = int_t(ty0)-iy0, y1 = int_t(ty1)-iy0;
                data(y0-_-y1,x0-_-x1) = std::move(subcut);
            } else {
                // Source fully covered
                int_t y0 = iy0, y1 = iy1, x0 = ix0, x1 = ix1;
                iimg.read_subset(data, y0-_-y1, x0-_-x1);
            }

            if (img.data.empty()) {
                // First image on which this source is found
                img.data = std::move(data);

                // Initialize WCS
                img.hdr = astro::filter_wcs(hdr);

                // Save precise center
                px = hs+1+(dxc-xc);
                py = hs+1+(dyc-yc);
            } else {
                // This source was found on another image, combine the two
                vec1u idb = where(!is_finite(img.data));
                img.data[idb] = data[idb];
            }
        }

        if (img.data.empty()) {
            // Source not covered, just use placeholder data
            img.data = replicate(dnan, 2*hs+1, 2*hs+1);
            img.hdr = astro::filter_wcs(first_hdr);
        }

        // Set cutout WCS
        if (!fits::setkey(img.hdr, "NAXIS1", img.data.dims[1]) ||
            !fits::setkey(img.hdr, "NAXIS2", img.data.dims[0]) ||
            !fits::setkey(img.hdr, "CRPIX1", px) || !fits::setkey(img.hdr, "CRPIX2", py) ||
            !fits::setkey(img.hdr, "CRVAL1", ra) || !fits::setkey(img.hdr, "CRVAL2", dec)) {
            write_warning("could not set WCS information (CRPIX1, CRPIX2, CRVAL1, CRVAL2)");
            write_warning("WCS for the cutout will be wrong");
            write_warning("parsing '"+img.source.short_name+"'");
        }

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
