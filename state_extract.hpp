void get_background_apertures(const vec2b& mask, double radius, vec1u& y, vec1u& x) {
    y.clear();
    x.clear();

    // First segment the mask
    uint_t nseg = 0;
    vec2u seg = segment(mask, nseg);

    // Add segments at the edges
    seg(0,_) = seg(seg.dims[0]-1,_) = seg(_,0) = seg(_,seg.dims[1]-1) = nseg+1;

    // Build distance from segments
    vec2d dist;
    vec2u dist_id;
    segment_distance(seg, dist, dist_id);

    // Now add circles on the largest distance point repeatedly until there is no
    // space left
    double radius2 = sqr(radius);
    uint_t mid = max_id(dist);
    while (dist[mid] >= radius2) {
        vec1u muid = mult_ids(dist, mid);
        y.push_back(muid[0]);
        x.push_back(muid[1]);
        dist[mid] = 0;
        dist_id[mid] = mid;

        std::vector<uint_t> oy = {muid[0]}, ox = {muid[1]};
        while (!ox.empty()) {
            auto toy = oy, tox = ox;
            oy.clear(); ox.clear();

            while (!tox.empty()) {
                uint_t ty = toy.back(); toy.pop_back();
                uint_t tx = tox.back(); tox.pop_back();

                auto check_add = [&](uint_t tty, uint_t ttx) {
                    if (dist.safe(tty,ttx) == 0) return;

                    if (dist_id.safe(ty,tx) == mid) {
                        double nd = sqr(double(tty) - y.back())
                                  + sqr(double(ttx) - x.back());

                        if (nd < radius2) {
                            // Still within the circle, move on
                            dist.safe(tty,ttx) = 0;
                            dist_id.safe(tty,ttx) = mid;
                            oy.push_back(tty);
                            ox.push_back(ttx);
                        } else {
                            // Setup new origin and continue below
                            dist_id.safe(ty,tx) = flat_id(dist, ty, tx);
                        }
                    }

                    if (dist_id.safe(ty,tx) != mid) {
                        double x0 = dist_id.safe(ty,tx) % dist_id.dims[1];
                        double y0 = dist_id.safe(ty,tx) / dist_id.dims[1];
                        double nd = sqr(double(tty) - y0) + sqr(double(ttx) - x0);
                        if (dist.safe(tty,ttx) > nd) {
                            dist.safe(tty,ttx) = nd;
                            dist_id.safe(tty,ttx) = dist_id.safe(ty,tx);
                            oy.push_back(tty);
                            ox.push_back(ttx);
                        }
                    }
                };

                if (ty != 0)              check_add(ty-1,tx);
                if (ty != dist.dims[0]-1) check_add(ty+1,tx);
                if (tx != 0)              check_add(ty,tx-1);
                if (tx != dist.dims[1]-1) check_add(ty,tx+1);
            }
        }

        mid = max_id(dist);
    }
}

void state_t::extract_fluxes() {
    if (!good) return;

    for (auto& img : images) {
        if (!img.covered || img.source.flux_method != "aper") continue;

        auto& rdata = img.regridded_data;
        auto& rpsf = img.regridded_psf;

        if (opts.verbose) note("extracting ", img.source.short_name);

        if (rdata.empty()) {
            // Regrid, if not done in the previous stages already
            if (opts.verbose) note("  regridding...");
            regrid(img);
        }

        if (!img.source.noconvolve && !img.convolved) {
            if (opts.verbose) note("  convolving for homogenization...");
            convolve_to(img, worst_seeing);
        }

        if (opts.verbose) note("  creating background apertures (", img.aperture, " arcsec)");

        // Create background regions
        double aper_radius = img.aperture/det_aspix/2.0;
        double threshold = opts.bg_threshold;
        vec1u bx, by;

        // Make sure the distance between each region is at least the aperture
        // radius, and if the image has been convolved, make sure it is also larger
        // than the convolution scale so that noise is not correlated.
        double min_dist = aper_radius;
        if (!img.source.noconvolve) {
            min_dist = sqrt(sqr(min_dist) + sqr(img.conv_radius));
        }

        // Flag center of image where the source is, just in case it didn't make it into
        // the detection image
        vec2b center = circular_mask(det_img.dims, max(det_seeing, img.seeing)/2.0/det_aspix,
            det_img.dims[0]/2, det_img.dims[1]/2) > 0.2;

        if (opts.no_neighbor_mask_background || img.source.no_neighbor_mask_background) {
            get_background_apertures(center || !is_finite(rdata), min_dist, by, bx);
        } else {
            segment_output sdo;
            segment_params sdp;
            sdp.min_area = 4;
            vec2b exclude = segment(det_img > threshold, sdo, sdp) > 0;
            get_background_apertures(exclude || center || !is_finite(rdata), min_dist, by, bx);

            if (by.size() < opts.min_bg_aper) {
                uint_t ntry = 0;
                while (by.size() < opts.min_bg_aper && ntry < 4) {
                    ++ntry;
                    write_warning("could place only ", by.size()," background apertures below ", threshold, " sigma");
                    threshold += 0.5;
                    exclude = segment(det_img > threshold, sdo, sdp) > 0;
                    write_warning("trying again with a threshold of ", threshold, " sigma");
                    get_background_apertures(exclude || center || !is_finite(rdata), min_dist, by, bx);
                }

                if (by.size() < opts.min_bg_aper) {
                    // Try without masking neighbors, just the central source
                    write_warning("could place only ", by.size()," background apertures below ", threshold, " sigma");
                    write_warning("trying again without masking neighbors");
                    get_background_apertures(center || !is_finite(rdata), min_dist, by, bx);
                }
            }
        }

        if (by.size() < opts.min_bg_aper) {
            write_error("could not place enough background apertures (needed ",
                opts.min_bg_aper, ", got ", by.size(), ")");
            write_error("please check the detection image");
            good = false;
            return;
        }

        // Write regions
        if (!opts.nocuts) {
            write_ds9_region(outdir+img.source.short_name+"_bg.reg",
                det_wcs, by, bx, aper_radius);
        }

        if (opts.verbose) note("  placed ", by.size(), " background apertures");

        // Create aperture mask
        vec2d aper_mask; {
            uint_t y0 = det_img.dims[0]/2;
            uint_t x0 = det_img.dims[1]/2;
            aper_mask = circular_mask(det_img.dims, aper_radius, y0, x0);

            // Write region
            if (!opts.nocuts) {
                write_ds9_region(outdir+img.source.short_name+"_aper.reg",
                    det_wcs, vec1u({y0}), vec1u({x0}), aper_radius);
            }
        }

        if (opts.verbose) {
            note("extracting aperture flux...");
        }

        // Compute centroid
        double y0 = 0, x0 = 0, taper = 0;
        for (uint_t iy : range(aper_mask.dims[0]))
        for (uint_t ix : range(aper_mask.dims[1])) {
            y0 += iy*aper_mask.safe(iy,ix);
            x0 += ix*aper_mask.safe(iy,ix);
            taper += aper_mask.safe(iy,ix);
        }

        y0 /= taper;
        x0 /= taper;

        // Find good pixels
        img.mask = vec2d(is_finite(rdata));
        rdata[where(img.mask < 0.5)] = 0;

        img.flux_bg = replicate(dnan, bx.size());
        vec1d bg_area(bx.size());
        for (uint_t i : range(bx)) {
            // Make background aperture mask
            vec2d bg_mask = circular_mask(aper_mask.dims, aper_radius, by[i], bx[i]);
            bg_area[i] = total(bg_mask);
            // Compute fractional area covered inside the mask
            double covcor = bg_area[i]/total(bg_mask*img.mask);
            // Sum pixels and correct for coverage
            img.flux_bg[i] = total(rdata*bg_mask*img.mask)*covcor;
        }

        // Clip strong background outliers
        vec1u idgbg = where(sigma_clip(img.flux_bg, 10.0));
        if (idgbg.size() < 3) {
            // Still not enough, just use everything
            idgbg = uindgen(img.flux_bg.size());
        }

        img.num_bg = idgbg.size();

        // Compute and subtract background level
        img.background = median(img.flux_bg[idgbg]/bg_area[idgbg]);
        rdata -= img.background;

        // Compute aperture flux
        img.flux = total(rdata*aper_mask*img.mask);

        // Compute uncertainty from background apertures
        // (the second term accounts for the uncertainty on the background level)
        img.flux_err = stddev(img.flux_bg[idgbg])*sqrt(1.0 + 1.0/idgbg.size());

        // Compute aperture correction from this image's PSF (assumes point source!)
        vec2d model = make_point_source(rpsf, rdata.dims, y0, x0);
        // vec2d model = translate(rpsf, y0 - rpsf.dims[0]/2, x0 - rpsf.dims[1]/2);
        img.apcor = 1.0/total(model*aper_mask*img.mask);

        // Apply aperture correction
        img.flux *= img.apcor;
        img.flux_err *= img.apcor;

        // Write regridded, convolved and background subtracted image
        if (!opts.nocuts) {
            rdata[where(img.mask < 0.5)] = dnan;
            fits::write(outdir+img.source.short_name+"_regrid.fits",
                rdata, img.regridded_hdr
            );
        }
    }
}
