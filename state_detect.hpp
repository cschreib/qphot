void get_background(const vec2d& data, double& bg, double& rms) {
    vec1u idg;
    idg = where(is_finite(data));
    bg = median(data[idg]);
    rms = 1.48*mad(data[idg]);

    idg = where(is_finite(data) && data - bg < 10*rms);
    bg = median(data[idg]);
    rms = 1.48*mad(data[idg]);
}

void state_t::build_detection_image() {
    if (!good) return;

    if (opts.verbose) {
        note("building detection image...");
    }

    double largest_size = 0;
    double ra0 = 0, dec0 = 0;
    uint_t ncovered = 0;
    bool has_detect = false;

    for (auto& img : images) {
        auto& idims = img.data.dims;

        // Apply zero point if needed
        if (is_finite(img.source.zero_point) && img.source.zero_point != 23.9) {
            img.data *= e10(0.4*(23.9 - img.source.zero_point));
        }

        // Mask bad pixels if asked

        // From DS9 region files
        std::string bad_file = outdir+img.source.short_name+"-bad.reg";
        if (!file::exists(bad_file)) {
            bad_file = outdir+img.source.short_name+"_bad.reg";
        }
        if (!file::exists(bad_file)) {
            bad_file = img.source.short_name+"-bad.reg";
        }
        if (!file::exists(bad_file)) {
            bad_file = img.source.short_name+"_bad.reg";
        }
        if (!file::exists(bad_file)) {
            bad_file = "";
        }

        if (!bad_file.empty()) {
            vec2d bad;
            try {
                read_ds9_region_circles_physical(bad_file, img.wcs, bad);
            } catch (std::runtime_error& e) {
                write_error(e.what());
                good = false;
                return;
            }

            vec2d bad_mask(idims);
            for (uint_t i : range(bad.dims[0])) {
                bad_mask += circular_mask(idims, bad(i,2), bad(i,1), bad(i,0));
            }

            img.data[where(bad_mask > 0.5)] = dnan;
        }

        // From FITS mask
        bad_file = outdir+img.source.short_name+"-bad.fits";
        if (!file::exists(bad_file)) {
            bad_file = outdir+img.source.short_name+"_bad.fits";
        }
        if (!file::exists(bad_file)) {
            bad_file = img.source.short_name+"-bad.fits";
        }
        if (!file::exists(bad_file)) {
            bad_file = img.source.short_name+"_bad.fits";
        }
        if (!file::exists(bad_file)) {
            bad_file = "";
        }

        if (!bad_file.empty()) {
            vec2b bad;
            fits::read(bad_file, bad);
            if (bad.dims == img.data.dims) {
                img.data[where(bad)] = dnan;
            }
        }

        // Check the content of the image
        // Mask portions of image only containing 0s (only if larger than 4 adjacent pixels)
        segment_output sdo;
        vec2u segz = segment(img.data == 0, sdo);
        foreach_segment(segz, sdo.origin, [&](uint_t is, vec1u ids) {
            if (sdo.area[is] > 4) {
                img.data[ids] = dnan;
            }
        });

        // Check fraction of non finite pixels
        double fnfin = fraction_of(!is_finite(img.data));
        if (fnfin > 0.4) {
            write_warning("ignoring image '", img.source.short_name, "' because it contains ",
                round(fnfin*100), "% of invalid values");
            img.data[_] = dnan;
            img.covered = false;
        }

        if (img.covered) {
            if (img.source.clean_neighbors) {
                has_clean = true;
            }
            if (!img.source.nodetect) {
                has_detect = true;
            }
            if (!img.source.noconvolve) {
                has_homogenize = true;
            }
            if (img.source.hri) {
                hri.push_back(&img);
            }

            if (img.source.aspix < smallest_aspix) {
                smallest_aspix = img.source.aspix;
            }

            double size = img.source.aspix*max(idims[0], idims[1]);
            if (size > largest_size) {
                largest_size = size;
            }

            double ra, dec;
            astro::xy2ad(img.wcs, idims[1]/2.0 + 1.0, idims[0]/2.0 + 1.0, ra, dec);
            ra0 += ra;
            dec0 += dec;

            if (!img.source.noconvolve) {
                if (img.seeing > worst_seeing) {
                    worst_seeing = img.seeing;
                }
            }

            if (!img.source.nodetect) {
                if (img.seeing > det_seeing) {
                    det_seeing = img.seeing;
                }
            }

            ++ncovered;
        }

        img.data_hash = hash(img.hdr, img.data, img.psf);
    }

    if (ncovered == 0) {
        write_error("no image provided valid data, cannot extract fluxes");
        good = false;
        return;
    }
    if (!has_detect) {
        write_error("no image provided valid data for detection, cannot extract fluxes");
        good = false;
        return;
    }
    if (opts.clean_model == "hri") {
        if (hri.size() > 1) {
            error("when opts.clean_model=hri, only one image can be flagged as hri=1 (got ", hri.size(), ")");
            good = false;
            return;
        } else if (hri.empty()) {
            error("no image found with hri=1, need one when opts.clean_model=hri");
            good = false;
            return;
        }
    }
    if (opts.clean_model == "nearest_hri") {
        if (hri.empty()) {
            error("no image found with hri=1, need at least one when opts.clean_model=nearest_hri");
            good = false;
            return;
        }
    }

    ra0 /= ncovered;
    dec0 /= ncovered;

    // Set aperture size if not defined
    if (!is_finite(aperture)) {
        aperture = worst_seeing;
    }

    for (auto& img : images) {
        if (!is_finite(img.aperture)) {
            img.aperture = aperture;
        }
    }

    // Set search radius if not defined
    if (!is_finite(search_radius)) {
        search_radius = det_seeing;
    }

    if (opts.verbose) {
        note("worst seeing for detection is ", det_seeing, "\"");
        note("worst seeing for homogenization is ", worst_seeing, "\"");
        note("using a default aperture of ", aperture, "\"");
    }

    // Check if we can reuse a cached detection image
    bool rebuild_detection = true;
    std::string det_outfile = outdir+"detection_image.fits";
    {
        vec1s hashes;
        for (auto& img : images) {
            if (!img.source.nodetect) {
                hashes.push_back(hash(img.data_hash, img.covered));
            }
        }

        det_cache_hash = hash(det_seeing, hashes);
    }

    if (opts.reuse && file::exists(det_outfile)) {
        fits::input_image iimg(det_outfile);

        std::string thash;
        if (iimg.read_keyword("CACHEID", thash)) {
            if (thash == det_cache_hash) {
                // We're good, read the cached image
                if (opts.verbose) {
                    note("reading detection image from ", det_outfile, " ...");
                }

                iimg.reach_hdu(1);
                iimg.read(det_img);
                det_hdr = iimg.read_header();
                det_wcs = astro::wcs(det_hdr);
                get_pixel_size(det_wcs, det_aspix);

                iimg.reach_hdu(2);
                iimg.read(det_psf);
                det_psf_hdr = iimg.read_header();
                det_psf_wcs = astro::wcs(det_psf_hdr);

                rebuild_detection = false;
            } else {
                if (opts.verbose) note("found a detection image but its CACHEID is different");
            }
        } else {
            if (opts.verbose) note("found a detection image but it does not have the CACHEID keyword");
        }
    }

    if (rebuild_detection) {
        if (opts.verbose) {
            note("defining stack WCS...");
        }

        // Make the combined image WCS
        {
            make_wcs_header_params p;
            // Get pixel scale from most refined image
            p.pixel_scale = det_aspix = smallest_aspix;
            p.sky_ref_ra = ra0;
            p.sky_ref_dec = dec0;
            // Get image size from largest image
            uint_t isize = ceil(largest_size/p.pixel_scale);
            // Make sure the size is odd
            if (isize % 2 == 0) ++isize;
            p.dims_x = p.dims_y = isize;
            p.pixel_ref_x = p.pixel_ref_y = p.dims_x/2.0 + 1.0;

            if (!make_wcs_header(p, det_hdr)) {
                write_error("could not make detection header WCS");
                good = false;
                return;
            }

            det_img.resize(p.dims_y, p.dims_x);

            // The PSF will be twice larger
            isize = 2*isize+1;
            p.dims_x = p.dims_y = isize;
            p.pixel_ref_x = p.pixel_ref_y = p.dims_x/2.0 + 1.0;

            if (!make_wcs_header(p, det_psf_hdr)) {
                write_error("could not make detection PSF header WCS");
                good = false;
                return;
            }

            det_psf.resize(p.dims_y, p.dims_x);
        }

        det_wcs = astro::wcs(det_hdr);
        det_psf_wcs = astro::wcs(det_psf_hdr);
        fits::setkey(det_hdr, "SEEING", det_seeing, "[arcsec]");
        fits::setkey(det_hdr, "CACHEID", "'"+det_cache_hash+"'");

        if (opts.verbose) {
            note("making detection image...");
        }

        // Regrid each image, convolve and co-add it
        vec3d cube;
        vec3d cube_psf;
        vec1d im_rms;
        for (auto& img : images) {
            if (!img.covered) {
                img.regridded_data = replicate(dnan, det_img.dims);
                continue;
            }

            if (img.source.nodetect) continue;

            if (opts.verbose) {
                note("handling ", img.source.short_name);
                note("  regridding...");
            }

            // Regrid
            regrid(img);

            // Make a copy, because we won't keep it
            vec2d tdata = img.regridded_data;
            vec2d tpsf = img.regridded_psf;

            // Convolve
            if (img.seeing < det_seeing) {
                if (opts.verbose) {
                    note("  convolving for detection...");
                }

                double conv_radius = sqrt(sqr(det_seeing) - sqr(img.seeing))/2.355/det_aspix;

                // Build kernel
                double ky0 = tdata.dims[0]/2;
                double kx0 = tdata.dims[1]/2;
                vec2d kernel = gaussian_profile(tdata.dims, conv_radius, ky0, kx0);

                // Convolve image
                vec1u idb = where(!is_finite(tdata));
                tdata[idb] = 0;
                tdata = convolve2d(tdata, kernel);

                idb = where(!is_finite(tpsf));
                tpsf[idb] = 0;
                tpsf = convolve2d(tpsf, kernel);
            }

            // Get RMS and background level
            if (opts.verbose) {
                note("  determining RMS and background...");
            }

            double med, rms;
            get_background(tdata, med, rms);
            im_rms.push_back(rms);

            // Coadd
            if (opts.verbose) {
                note("  coadding...");
            }

            append<2>(cube, reform(tdata - med, tdata.dims, 1));
            append<2>(cube_psf, reform(tpsf, tpsf.dims, 1));
        }

        // Build stacked image and add it to the cube
        vec2d stack; {
            vec1d wei = 1/sqr(im_rms);
            double twei = total(wei);

            stack = reduce(2, cube, [&](vec1d v) {
                return total(v*wei)/twei;
            });

            double md;
            im_rms.push_back(0.0);
            get_background(stack, md, im_rms.back());
            stack -= md;
            append<2>(cube, reform(stack, stack.dims, 1));
        }

        // Find maximum SNR
        det_img = reduce(2, cube, [&](vec1d v) {
            return max(v/im_rms);
        });
        det_psf = partial_mean(2, cube_psf);

        // Write detection image
        if (!opts.nocuts) {
            fits::output_image oimg(det_outfile);
            oimg.reach_hdu(1);
            oimg.write(det_img);
            oimg.write_header(det_hdr);
            oimg.reach_hdu(2);
            oimg.write(det_psf);
            oimg.write_header(det_psf_hdr);
        }
    }
}
