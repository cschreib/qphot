void state_t::apply_filters() {
    if (!good || !has_filters) return;

    if (opts.verbose) {
        note("applying filters...");
    }

    for (auto& img : images) {
        if (!img.covered) continue;

        if (img.source.background_box > 0) {
            if (opts.verbose) note("  applying background filter to ", img.source.short_name, "...");

            auto& rdata = img.regridded_data;
            if (rdata.empty()) {
                // Regrid, if not done at the detection stage already
                if (opts.verbose) note("  regridding...");
                regrid(img);
            }

            if (!img.source.noconvolve && !img.convolved) {
                if (opts.verbose) note("  convolving for homogenization...");
                convolve_to(img, worst_seeing);
            }

            uint_t dsize = img.source.background_box*0.5/det_aspix;
            vec2b mask = mask_inflate(seg > 0, 2);

            // Compute background in a moving box using median on the image, masking sources
            vec2d background = replicate(dnan, rdata.dims);
            for (uint_t iy : range(rdata.dims[0]))
            for (uint_t ix : range(rdata.dims[1])) {
                uint_t y0 = (iy > dsize ? iy - dsize : 0);
                uint_t x0 = (ix > dsize ? ix - dsize : 0);
                uint_t y1 = (iy < rdata.dims[0]-dsize ? iy + dsize : rdata.dims[0]-1);
                uint_t x1 = (ix < rdata.dims[1]-dsize ? ix + dsize : rdata.dims[1]-1);

                vec2d nb = rdata(y0-_-y1,x0-_-x1);
                vec2b sb = mask(y0-_-y1,x0-_-x1);
                vec1u idg = where(!sb);
                if (idg.size() > 4) {
                    background(iy,ix) = median(nb[idg]);
                }
            }

            // Fill in the gaps where background could not be determined using nearby values
            // Use zero as a fallback, i.e. assume background was ok
            vec2d old_background = background;
            for (uint_t iy : range(rdata.dims[0]))
            for (uint_t ix : range(rdata.dims[1])) {
                if (is_finite(old_background(iy,ix))) continue;

                uint_t y0 = (iy > dsize ? iy - dsize : 0);
                uint_t x0 = (ix > dsize ? ix - dsize : 0);
                uint_t y1 = (iy < rdata.dims[0]-dsize ? iy + dsize : rdata.dims[0]-1);
                uint_t x1 = (ix < rdata.dims[1]-dsize ? ix + dsize : rdata.dims[1]-1);

                vec2d nb = old_background(y0-_-y1,x0-_-x1);
                vec1u idg = where(is_finite(nb));
                if (idg.size() > 0) {
                    background(iy,ix) = mean(nb[idg]);
                } else {
                    background(iy,ix) = 0;
                }
            }

            // Subtract background from image
            rdata -= background;
        }
    }
}
