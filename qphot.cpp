#include <phypp.hpp>
#include "log.hpp"
#include "image.hpp"
#include "ds9.hpp"

void get_background(const vec2d& data, double& bg, double& rms) {
    vec1u idg;
    idg = where(is_finite(data));
    bg = median(data[idg]);
    rms = 1.48*mad(data[idg]);

    idg = where(is_finite(data) && data - bg < 10*rms);
    bg = median(data[idg]);
    rms = 1.48*mad(data[idg]);
}

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

int phypp_main(int argc, char* argv[]) {
    std::string outdir = "./";
    bool verbose = false;
    double aperture = dnan; // [arcsec]
    std::string exclude;
    std::string include;
    double det_seeing = 0;
    uint_t det_min_area = 30;
    std::string error_file;
    uint_t min_bg_aper = 10;
    bool no_neighbor_mask_background = false;
    bool reuse = false;
    std::string clean_model = "clean";
    double clean_ratio = 0.1;
    double clean_threshold = 3.0;
    double clean_conv_seeing = 0.4;
    bool save_models = false;
    bool debug_clean = false;
    bool clean_group = false;
    bool model_dance = false;
    double dance_step = 0.1;
    double dance_max = 0.5;
    double dance_chi2_range = 6.0;

    read_args(argc-1, argv+1, arg_list(
        aperture, outdir, verbose, exclude, include, det_seeing, error_file, min_bg_aper,
        clean_model, reuse, clean_ratio, clean_threshold, clean_conv_seeing, clean_group,
        model_dance, dance_step, dance_max, dance_chi2_range, det_min_area, debug_clean,
        no_neighbor_mask_background, save_models
    ));

    if (!error_file.empty()) {
        errfile.open(error_file, std::ios::app);
    }

    if (clean_model != "psf" && clean_model != "hri" && clean_model != "nearest_hri" &&
        clean_model != "clean") {
        error("unknown clean model '", clean_model, "'");
        note("allowed values: psf, hri, nearest_hri, clean");
        return 1;
    }

    outdir = file::directorize(outdir);
    file::mkdir(outdir);

    // Read filter data base
    auto fdb = read_filter_db(data_dir+"fits/filter-db/db.dat");
    auto fmap = read_filter_map(data_dir+"fits/filter-db/fast.map");

    // Get files and their properties
    vec<1,image_t> images;
    if (!read_image_list(argv[1], images)) {
        return 1;
    }

    // Remove bands we don't want
    if (!exclude.empty() || !include.empty()) {
        vec1u idb;
        for (uint_t i : range(images)) {
            if ((!exclude.empty() && regex_match(images[i].filename, exclude)) ||
                (!include.empty() && !regex_match(images[i].filename, include))) {
                idb.push_back(i);
            }
        }

        inplace_remove(images, idb);
    }

    if (verbose) {
        note("reading input images...");
    }

    double smallest_aspix = dinf;
    double largest_size = 0;
    double worst_seeing = 0;
    double ra0 = 0, dec0 = 0;
    uint_t ncovered = 0;
    uint_t nhri = 0;
    image_t* hri = nullptr;
    bool has_clean = false;
    bool has_detect = false;
    bool has_homogenize = false;

    for (auto& img : images) {
        img.input_image = std::unique_ptr<fits::input_image>(
            new fits::input_image(img.filename)
        );

        auto& iimg = *img.input_image;

        if (!iimg.is_image()) {
            write_error("'", img.filename, "' is not a FITS image");
            return 1;
        }

        img.hdr = iimg.read_header();
        img.wcs = astro::wcs(img.hdr);
        if (!astro::get_pixel_size(img.wcs, img.aspix)) {
            write_error("could not find pixel size in '", img.filename, "'");
            return 1;
        }

        filter_t fil;
        if (!get_filter(fdb, img.band, fil)) {
            write_warning("unknown filter '", img.band, "'");
        } else {
            uint_t eid = npos;
            get_filter_id(fmap, img.band, eid);

            img.lambda = fil.rlam;
            img.eazy_band = eid;
        }

        iimg.read(img.data);

        if (is_finite(img.zero_point)) {
            img.data *= e10(0.4*(23.9 - img.zero_point));
        }

        if (!img.psffile.empty()) {
            // Read the PSF
            fits::read(img.psffile, img.psf);
            if (!is_finite(img.seeing)) {
                img.seeing = get_seeing(img.psf, img.psf.dims[0]/2, img.psf.dims[1]/2);
            }

            // Adapt the PSF image to the dimensions of the analyzed image
            // and make sure it is centered
            vec1i idm = mult_ids(img.psf, max_id(img.psf));
            img.psf = recenter(img.psf, idm[0], idm[1], img.data.dims);
        } else {
            // Create a PSF if it is not provided
            img.psf = gaussian_profile(img.data.dims, img.seeing/2.355/img.aspix,
                img.data.dims[0]/2, img.data.dims[1]/2
            );
        }

        double fzero = fraction_of(img.data == 0);
        double fnfin = fraction_of(!is_finite(img.data));
        if (fzero > 0.4) {
            write_warning("ignoring image '", img.filename, "' because it contains ",
                round(fzero*100), "% of zero values");
            img.data[_] = dnan;
            img.covered = false;
        } else if (fnfin > 0.4) {
            write_warning("ignoring image '", img.filename, "' because it contains ",
                round(fnfin*100), "% of invalid values");
            img.data[_] = dnan;
            img.covered = false;
        }

        if (img.covered) {
            if (img.clean_neighbors) {
                has_clean = true;
            }
            if (!img.nodetect) {
                has_detect = true;
            }
            if (!img.noconvolve) {
                has_homogenize = true;
            }
            if (img.hri) {
                hri = &img;
                ++nhri;
            }
        }

        if (img.covered) {
            if (img.aspix < smallest_aspix) {
                smallest_aspix = img.aspix;
            }

            double size = img.aspix*max(img.data.dims[0], img.data.dims[1]);
            if (size > largest_size) {
                largest_size = size;
            }

            double ra, dec;
            astro::xy2ad(img.wcs,
                img.data.dims[1]/2.0 + 1.0, img.data.dims[0]/2.0 + 1.0,
                ra, dec
            );
            ra0 += ra;
            dec0 += dec;

            if (!img.noconvolve) {
                if (img.seeing > worst_seeing) {
                    worst_seeing = img.seeing;
                }
            }

            if (!img.nodetect) {
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
        return 1;
    }
    if (!has_detect) {
        write_error("no image provided valid data for detection, cannot extract fluxes");
        return 1;
    }
    if (clean_model == "hri") {
        if (nhri > 1) {
            error("when clean_model=hri, only one image can be flagged as hri=1 (got ", nhri, ")");
            return 1;
        } else if (nhri == 0) {
            error("no image found with hri=1, need one when clean_model=hri");
            return 1;
        }
    }
    if (clean_model == "nearest_hri") {
        if (nhri == 0) {
            error("no image found with hri=1, need at least one when clean_model=nearest_hri");
            return 1;
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

    if (verbose) {
        note("worst seeing for detection is ", det_seeing, "\"");
        note("worst seeing for homogenization is ", worst_seeing, "\"");
        note("using a default aperture of ", aperture, "\"");
    }

    // Check if we can reuse a cached detection image
    bool rebuild_detection = true;
    std::string det_outfile = outdir+"detection_image.fits";
    std::string det_cache_hash; {
        vec1s hashes;
        for (auto& img : images) {
            if (!img.nodetect) {
                hashes.push_back(hash(img.data_hash, img.covered));
            }
        }

        det_cache_hash = hash(det_seeing, hashes);
    }

    vec2d det_img;
    vec2d det_psf;
    astro::wcs det_wcs;
    double det_aspix = 0;
    fits::header det_hdr;
    if (reuse && file::exists(det_outfile)) {
        fits::input_image iimg(det_outfile);

        std::string thash;
        if (iimg.read_keyword("CACHEID", thash)) {
            if (thash == det_cache_hash) {
                // We're good, read the cached image
                if (verbose) {
                    note("reading detection image from ", det_outfile, " ...");
                }

                iimg.reach_hdu(1);
                iimg.read(det_img);
                det_hdr = iimg.read_header();
                det_wcs = astro::wcs(det_hdr);
                get_pixel_size(det_wcs, det_aspix);

                iimg.reach_hdu(2);
                iimg.read(det_psf);

                rebuild_detection = false;
            } else {
                note("found a detection image but its CACHEID is different");
            }
        } else {
            note("found a detection image but it does not have the CACHEID keyword");
        }
    }

    if (rebuild_detection) {
        if (verbose) {
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
                return 1;
            }

            det_img.resize(p.dims_y, p.dims_x);
            det_psf.resize(p.dims_y, p.dims_x);
        }

        det_wcs = astro::wcs(det_hdr);
        fits::setkey(det_hdr, "SEEING", det_seeing, "[arcsec]");
        fits::setkey(det_hdr, "CACHEID", "'"+det_cache_hash+"'");

        if (verbose) {
            note("making detection image...");
        }

        // Regrid each image, convolve and co-add it
        vec3d cube;
        vec3d cube_psf;
        vec1d im_rms;
        for (uint_t i : range(images)) {
            auto& img = images[i];
            auto& rdata = img.regridded_data;

            if (!img.covered) {
                rdata = replicate(dnan, det_img.dims);
                continue;
            }

            if (img.nodetect) continue;

            if (verbose) {
                note("handling ", img.filename);
                note("  regridding...");
            }

            // Regrid
            img.regrid_data(det_wcs, det_hdr, det_aspix, reuse, det_cache_hash, outdir);

            // Make a copy, because we won't keep it
            vec2d tdata = rdata;
            vec2d tpsf = img.regridded_psf;

            // Convolve
            if (img.seeing < det_seeing) {
                if (verbose) {
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
            if (verbose) {
                note("  determining RMS and background...");
            }

            double med, rms;
            get_background(tdata, med, rms);
            im_rms.push_back(rms);

            // Coadd
            if (verbose) {
                note("  coadding...");
            }

            append<2>(cube, reform(tdata - med, tdata.dims, 1));
            append<2>(cube_psf, reform(tpsf, tdata.dims, 1));
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
        fits::output_image oimg(det_outfile);
        oimg.reach_hdu(1);
        oimg.write(det_img);
        oimg.write_header(det_hdr);
        oimg.reach_hdu(2);
        oimg.write(det_psf);
    }

    if (has_clean) {
        if (verbose) {
            note("cleaning images...");
        }

        if (verbose) {
            note("cleaning detection image...");
        }

        // vec2d det_psf;

        // if () {
        //     // Only one image used for detection, use it's PSF directly
        //     det_psf = dimg.regridded_psf;
        // } else {
        //     // Build an effective PSF from the chosen detection seeing
        //     det_psf = gaussian_profile(det_img.dims,
        //         det_seeing/2.355/det_aspix, det_img.dims[0]/2, det_img.dims[1]/2
        //     );

        //     det_psf /= max(det_psf);
        // }

        vec2d model(det_img.dims);
        vec2d det_img_clean; {
            vec2d tdimg = det_img;
            uint_t mid = max_id(det_img);
            uint_t iter = 0;
            while (tdimg[mid] > clean_threshold) {
                vec1i ids = mult_ids(tdimg, mid);
                vec2d tpsf = translate_integer(det_psf,
                    ids[0] - int_t(det_psf.dims[0])/2, ids[1] - int_t(det_psf.dims[1])/2
                );

                if (debug_clean && iter % 20 == 0) {
                    print("clean ", iter, ": max ", tdimg[mid], " pos: ", ids);
                    if (iter % 100 == 0) {
                        fits::write("debug_clean.fits", tdimg);
                    }
                }

                double flx = tdimg[mid]*clean_ratio;
                tdimg -= tpsf*flx;
                model[mid] += flx;

                mid = max_id(tdimg);

                ++iter;
            }

            vec2d model_psf = gaussian_profile(det_img.dims,
                clean_conv_seeing/2.355/det_aspix, det_img.dims[0]/2, det_img.dims[1]/2
            );

            model_psf /= max(model_psf);
            det_img_clean = convolve2d(model, model_psf);

            fits::output_image oimg(outdir+"detection_image_cleaned.fits");
            oimg.reach_hdu(1);
            oimg.write(det_img_clean);
            oimg.write_header(det_hdr);
            oimg.write_keyword("CRATIO", clean_ratio);
            oimg.write_keyword("CTHRESH", clean_threshold);
            oimg.write_keyword("CSEING", clean_conv_seeing);

            oimg.reach_hdu(2);
            oimg.write(model);
            oimg.write_header(det_hdr);
        }

        if (verbose) {
            note("building segmentation map...");
        }

        segment_deblend_params sdp;
        sdp.detect_threshold = 0.5;
        sdp.deblend_threshold = det_seeing/det_aspix;
        sdp.min_area = det_min_area;
        segment_deblend_output sdo;
        vec2u seg = segment_deblend(det_img_clean, sdo, sdp);

        uint_t isource = npos; {
            // Search among the detected segment for a source close to the center
            double bestd = finf;
            for (uint_t i : range(sdo.id)) {
                double d = sqr(sdo.bx[i] - int_t(seg.dims[1]/2)) +
                           sqr(sdo.by[i] - int_t(seg.dims[0]/2));
                if (d < bestd) {
                    bestd = d;
                    isource = i;
                }
            }

            if (bestd > sqr(det_seeing/2.355/det_aspix)) {
                write_warning("the target source was not found in the detection image");
                write_warning("(the closest source was found at ", sqrt(bestd)*det_aspix, " arcsec)");
                write_warning("it will be added as a central point source for the fitting stage");

                isource = sdo.bx.size();
                sdo.id.push_back(max(sdo.id)+1);
                sdo.px.push_back(det_img.dims[1]/2);
                sdo.py.push_back(det_img.dims[0]/2);
                sdo.bx.push_back(sdo.px.back());
                sdo.by.push_back(sdo.py.back());
                sdo.area.push_back(1);
                sdo.origin.push_back(flat_id(seg, sdo.py.back(), sdo.px.back()));
                sdo.flux.push_back(det_img(sdo.py.back(), sdo.px.back()));

                seg(sdo.py.back(), sdo.px.back()) = sdo.id.back();
                model(sdo.py.back(), sdo.px.back()) = 1.0;
            }
        }

        if (verbose) {
            note("found ", sdo.id.size(), " segments in the detection image");
        }

        // Save segmentation map
        fits::write(outdir+"segmentation.fits", seg, det_hdr);
        fits::write_table(outdir+"segmentation_catalog.fits", ftable(
            sdo.id, sdo.px, sdo.py, sdo.bx, sdo.by, sdo.area, sdo.origin, sdo.flux
        ));

        // Clean homogenized images
        if (has_homogenize && clean_group) {
            if (verbose) note("cleaning homogenized images...");

            // Make sure the images are homogenized first
            vec2b fitmask = replicate(true, det_img.dims);
            uint_t nhomo = 0;
            for (auto& img : images) {
                if (!img.covered || !img.clean_neighbors || img.noconvolve) continue;

                if (verbose) note("  preparing ", img.filename);

                auto& rdata = img.regridded_data;
                if (rdata.empty()) {
                    // Regrid, if not done at the detection stage already
                    if (verbose) note("  regridding...");
                    img.regrid_data(det_wcs, det_hdr, det_aspix, reuse, det_cache_hash, outdir);
                }

                if (!img.noconvolve && !img.convolved) {
                    if (verbose) note("  convolving for homogenization...");
                    img.convolve_to(worst_seeing);
                }

                fitmask = fitmask && is_finite(rdata);
                img.make_finite();
                ++nhomo;
            }

            // Build models
            if (verbose) note("  building models...");
            vec2d models(1 + sdo.id.size(), det_img.dims[0]*det_img.dims[1]);
            models(0,_) = 1.0; // background
            vec1u gmodel; gmodel.push_back(0);

            if (clean_model == "psf") {
                // Assume point sources
                for (uint_t s : range(sdo.id)) {
                    models(s+1,_) = flatten(translate(det_psf,
                        sdo.py[s] - int_t(det_img.dims[0]/2), sdo.px[s] - int_t(det_img.dims[1]/2)
                    ));
                }
            } else if (clean_model == "clean") {
                uint_t imod = 1;
                foreach_segment(seg, sdo.origin, [&](const vec1u& ids) {
                    // Get pixels of that segment
                    vec2d tmod(det_img.dims);
                    tmod[ids] = model[ids]/max(model[ids]);
                    // Convolve to the resolution of the fitted image
                    tmod = convolve2d(tmod, det_psf);
                    // Save model
                    models(imod,_) = flatten(tmod);
                    ++imod;
                });
            } else if (clean_model == "nearest_hri" || clean_model == "hri") {
                write_error("HRI fitting is not implemented");
                return 1;
            }

            // Inspect models and flag out the ones that are zero
            for (uint_t s : range(sdo.id)) {
                if (total(sqr(models(s+1,_))) > 0) {
                    gmodel.push_back(s+1);
                } else {
                    if (s == isource) {
                        write_error("fit model for the target only contains zero values");
                        return 1;
                    } else {
                        write_warning("model for source ", s, " is zero");
                        write_warning("this source will be ignored in the fit");
                    }
                }
            }

            vec1u idf = where(fitmask);
            vec2f fit_fluxes(nhomo,gmodel.size());
            vec1f fit_chi2 = replicate(finf, nhomo);
            vec1i fit_dx = replicate(0, nhomo);
            vec1i fit_dy = replicate(0, nhomo);
            vec1f err = replicate(1.0, idf.size());
            auto batch = linfit_pack_batch<vec1d>(err, models(gmodel,idf));
            if (!batch.fr.success) {
                write_error("could not perform the fit to remove the neighbors");
                return 1;
            }

            int_t dance_step_pix = max(1, ceil(dance_step/det_aspix));
            int_t dance_nstep = ceil(dance_max/det_aspix/dance_step_pix);
            auto dofit = [&](int_t dy, int_t dx) {
                uint_t ihomo = 0;
                for (auto& img : images) {
                    if (!img.covered || !img.clean_neighbors || img.noconvolve) continue;

                    img.fit_neighbors(batch, models, idf, gmodel, dance_chi2_range,
                        dance_step_pix*dy, dance_step_pix*dx
                    );

                    // See if this has improved
                    if (batch.fr.chi2 < fit_chi2[ihomo]) {
                        fit_chi2[ihomo] = batch.fr.chi2;
                        fit_fluxes(ihomo,_) = batch.fr.params;
                        fit_dx[ihomo] = dx;
                        fit_dy[ihomo] = dy;
                    }

                    ++ihomo;
                }
            };

            if (verbose) note("  fitting...");

            // Do the fit
            dofit(0, 0);

            if (model_dance && dance_nstep > 0) {
                if (verbose) {
                    note("  dancing (", 2*dance_nstep+1, "x", 2*dance_nstep+1, ", step ",
                        dance_step_pix*det_aspix, "\")...");
                }

                // Adjust positions
                vec2d tmodels = models;
                for (int_t dy = -dance_nstep; dy <= dance_nstep; ++dy)
                for (int_t dx = -dance_nstep; dx <= dance_nstep; ++dx) {
                    if (dx == 0 && dy == 0) {
                        // We already did that
                        continue;
                    }

                    dofit(dy, dx);
                }

                uint_t ihomo = 0;
                for (auto& img : images) {
                    if (!img.covered || !img.clean_neighbors || img.noconvolve) continue;

                    if (fit_dx[ihomo] != 0 || fit_dy[ihomo] != 0) {
                        // Shift the image for flux measurement later on
                        img.regridded_data = translate_integer(img.regridded_data,
                            -dance_step_pix*fit_dy[ihomo], -dance_step_pix*fit_dx[ihomo]
                        );
                    }

                    fits::setkey(img.regridded_hdr, "DANCEX", img.regridded_aspix*dance_step_pix*fit_dx[ihomo], "dance offset [arcsec]");
                    fits::setkey(img.regridded_hdr, "DANCEY", img.regridded_aspix*dance_step_pix*fit_dy[ihomo], "dance offset [arcsec]");
                    fits::setkey(img.regridded_hdr, "DANCES", img.regridded_aspix*dance_step_pix, "dance step [arcsec]");
                    fits::setkey(img.regridded_hdr, "DANCER", img.regridded_aspix*dance_nstep*dance_step_pix, "dance radius [arcsec]");

                    ++ihomo;
                }
            }

            if (verbose) note("  subtracting neighbors...");

            // Remove all the fitted segments from the image
            uint_t ihomo = 0;
            for (auto& img : images) {
                if (!img.covered || !img.clean_neighbors || img.noconvolve) continue;

                for (uint_t m : range(gmodel)) {
                    if (gmodel[m] == isource+1) {
                        // Keep our target there
                        continue;
                    }

                    if (!is_finite(fit_fluxes(ihomo,m))) {
                        write_warning("in the fit for image ", img.filename);
                        if (gmodel[m] == 0) {
                            write_warning("the fitted value of the background was invalid");
                            write_warning("something is very wrong!");
                        } else {
                            write_warning("the fitted value of source ", gmodel[m]-1, " was invalid");
                        }
                        continue;
                    }

                    vec2d tmod = reform(models(gmodel[m],_), img.regridded_data.dims);
                    img.regridded_data -= fit_fluxes(ihomo,m)*tmod;
                }

                img.restore_not_finite();
                img.cleaned = true;
                ++ihomo;
            }
        }

        for (auto& img : images) {
            if (!img.covered || !img.clean_neighbors || img.cleaned) continue;

            if (verbose) note("cleaning ", img.filename);

            auto& rdata = img.regridded_data;
            auto& rpsf = img.regridded_psf;

            if (rdata.empty()) {
                // Regrid, if not done at the detection stage already
                if (verbose) note("  regridding...");
                img.regrid_data(det_wcs, det_hdr, det_aspix, reuse, det_cache_hash, outdir);
            }

            if (!img.noconvolve && !img.convolved) {
                if (verbose) note("  convolving for homogenization...");
                img.convolve_to(worst_seeing);
            }

            if (verbose) note("  building models...");

            // Build models
            vec2d models(1 + sdo.id.size(), rdata.dims[0]*rdata.dims[1]);
            models(0,_) = 1.0; // background
            vec1u gmodel; gmodel.push_back(0);

            if (clean_model == "psf") {
                // Assume point sources
                for (uint_t s : range(sdo.id)) {
                    models(s+1,_) = flatten(translate(rpsf,
                        sdo.py[s] - int_t(rdata.dims[0]/2), sdo.px[s] - int_t(rdata.dims[1]/2)
                    ));
                }
            } else if (clean_model == "clean") {
                uint_t imod = 1;
                foreach_segment(seg, sdo.origin, [&](const vec1u& ids) {
                    // Get pixels of that segment
                    vec2d tmod(det_img.dims);
                    tmod[ids] = model[ids]/max(model[ids]);
                    // Convolve to the resolution of the fitted image
                    tmod = convolve2d(tmod, img.regridded_psf);
                    // Save model
                    models(imod,_) = flatten(tmod);
                    ++imod;
                });
            } else if (clean_model == "nearest_hri" || clean_model == "hri") {
                write_error("HRI fitting is not implemented");
                return 1;
            }

            // Inspect models and flag out the ones that are zero
            for (uint_t s : range(sdo.id)) {
                if (total(sqr(models(s+1,_))) > 0) {
                    gmodel.push_back(s+1);
                } else {
                    if (s == isource) {
                        write_error("fit model for the target only contains zero values");
                        return 1;
                    } else {
                        write_warning("model for source ", s, " is zero");
                        write_warning("this source will be ignored in the fit");
                    }
                }
            }

            if (verbose) note("  fitting...");

            // Do the fit
            img.make_finite();
            vec1u idf = img.idf;
            vec1d err = replicate(1.0, idf.size());
            vec1d fit_fluxes(gmodel.size());
            int_t fit_dx = 0, fit_dy = 0;
            float fit_chi2 = finf;
            auto batch = linfit_pack_batch<vec1d>(err, models(gmodel,idf));
            if (!batch.fr.success) {
                write_error("could not perform the fit to remove the neighbors");
                return 1;
            }

            int_t dance_step_pix = max(1, ceil(dance_step/img.regridded_aspix));
            int_t dance_nstep = ceil(dance_max/img.regridded_aspix/dance_step_pix);
            auto dofit = [&](int_t dy, int_t dx) {
                img.fit_neighbors(batch, models, idf, gmodel, dance_chi2_range,
                    dance_step_pix*dy, dance_step_pix*dx
                );

                if (batch.fr.success && batch.fr.chi2 < fit_chi2) {
                    fit_chi2 = batch.fr.chi2;
                    fit_fluxes = batch.fr.params;
                    fit_dy = dy;
                    fit_dx = dx;
                }
            };

            dofit(0, 0);

            // Adjust positions if asked
            if (model_dance && dance_nstep > 0) {
                if (verbose) {
                    note("  dancing (", 2*dance_nstep+1, "x", 2*dance_nstep+1, ", step ",
                        dance_step_pix*det_aspix, "\")...");
                }

                for (int_t dy = -dance_nstep; dy <= dance_nstep; ++dy)
                for (int_t dx = -dance_nstep; dx <= dance_nstep; ++dx) {
                    if (dx == 0 && dy == 0) {
                        // We already did that
                        continue;
                    }

                    dofit(dy, dx);
                }

                if (fit_dx != 0 || fit_dy != 0) {
                    // Shift the image for flux measurement later on
                    img.regridded_data = translate_integer(img.regridded_data,
                        -dance_step_pix*fit_dy, -dance_step_pix*fit_dx
                    );
                }

                fits::setkey(img.regridded_hdr, "DANCEX", img.regridded_aspix*dance_step_pix*fit_dx, " [arcsec]");
                fits::setkey(img.regridded_hdr, "DANCEY", img.regridded_aspix*dance_step_pix*fit_dy, " [arcsec]");
                fits::setkey(img.regridded_hdr, "DANCES", img.regridded_aspix*dance_step_pix, " [arcsec]");
                fits::setkey(img.regridded_hdr, "DANCER", img.regridded_aspix*dance_nstep*dance_step_pix, " [arcsec]");
            }

            if (verbose) note("  subtracting neighbors...");

            // Remove all the fitted segments from the image
            for (uint_t m : range(gmodel)) {
                if (gmodel[m] == isource+1) {
                    // Keep our target there
                    continue;
                }

                if (!is_finite(fit_fluxes[m])) {
                    if (gmodel[m] == 0) {
                        write_warning("the fitted value of the background was invalid");
                        write_warning("something is very wrong!");
                    } else {
                        write_warning("the fitted value of source ", gmodel[m]-1, " was invalid");
                    }
                    continue;
                }

                vec2d tmod = reform(models(gmodel[m],_), rdata.dims);
                rdata -= fit_fluxes[m]*tmod;
            }

            if (save_models) {
                vec2d tmodel(rdata.dims);
                for (uint_t m : range(gmodel)) {
                    if (!is_finite(fit_fluxes[m])) continue;

                    vec2d tmod = reform(models(gmodel[m],_), rdata.dims);
                    tmodel += fit_fluxes[m]*tmod;
                }

                fits::output_image oimg(outdir+img.base_filename+"_regrid_model.fits");
                oimg.write(tmodel);
                oimg.write_header(det_hdr);

                // Make catalog
                uint_t nsrc = count(gmodel != 0);
                vec1u id(nsrc);
                vec1d sx(nsrc), sy(nsrc), flux(nsrc);

                uint_t im = 0;
                for (uint_t m : range(gmodel)) {
                    if (gmodel[m] == 0) continue;

                    flux[im] = fit_fluxes[m];
                    id[im] = sdo.id[gmodel[m]-1];
                    sx[im] = sdo.px[gmodel[m]-1];
                    sy[im] = sdo.py[gmodel[m]-1];

                    ++im;
                }
                vec1u idsource = gmodel[where(gmodel != 0)];

                fits::output_table otbl(outdir+img.base_filename+"_regrid_model_catalog.fits");
                otbl.reach_hdu(1);
                otbl.write_columns("id", id, "x", sx, "y", sy, "flux", flux);
            }

            img.restore_not_finite();
            img.cleaned = true;
        }
    }

    for (auto& img : images) {
        if (!img.covered || img.flux_method != "aper") continue;

        auto& rdata = img.regridded_data;
        auto& rpsf = img.regridded_psf;

        if (verbose) note("extracting ", img.filename);

        if (rdata.empty()) {
            // Regrid, if not done in the previous stages already
            if (verbose) note("  regridding...");
            img.regrid_data(det_wcs, det_hdr, det_aspix, reuse, det_cache_hash, outdir);
        }

        if (!img.noconvolve && !img.convolved) {
            if (verbose) note("  convolving for homogenization...");
            img.convolve_to(worst_seeing);
        }

        if (verbose) note("  creating background apertures (", img.aperture, " arcsec)");

        // Create background regions
        double bg_radius = img.aperture/det_aspix/2.0;
        double threshold = 2.0;
        vec1u bx, by;

        // Make sure the distance between each region is at least the aperture
        // radius, and if the image has been convolved, make sure it is also larger
        // than the convolution scale so that noise is not correlated.
        double min_dist = bg_radius;
        if (!img.noconvolve) {
            min_dist = sqrt(sqr(min_dist) + sqr(img.conv_radius));
        }

        // Flag center of image where the source is, just in case it didn't make it into
        // the detection image
        vec2b center(det_img.dims);
        center(det_img.dims[0]/2, det_img.dims[1]/2) = true;

        if (no_neighbor_mask_background) {
            get_background_apertures(center || !is_finite(rdata), min_dist, by, bx);
        } else {
            get_background_apertures(det_img > threshold || center || !is_finite(rdata), min_dist, by, bx);

            if (by.size() < min_bg_aper) {
                uint_t ntry = 0;
                while (by.size() < min_bg_aper && ntry < 4) {
                    ++ntry;
                    write_warning("could place only ", by.size()," background apertures below ", threshold, " sigma");
                    threshold += 0.5;
                    write_warning("trying again with a threshold of ", threshold, " sigma");
                    get_background_apertures(det_img > threshold || center || !is_finite(rdata), min_dist, by, bx);
                }

                if (by.size() < min_bg_aper) {
                    write_warning("could place only ", by.size()," background apertures below ", threshold, " sigma");
                    threshold += 0.5;
                    write_warning("trying again without masking neighbors");

                    // Try without masking neighbors, just the central source
                    get_background_apertures(center || !is_finite(rdata), min_dist, by, bx);
                }
            }
        }

        if (by.size() < min_bg_aper) {
            write_error("could not place enough background apertures (needed ",
                min_bg_aper, ", got ", by.size(), ")");
            write_error("please check the detection image");
            return 1;
        }

        // Write regions
        write_ds9_region(outdir+img.base_filename+"_bg.reg",
            det_wcs, by, bx, bg_radius);

        if (verbose) note("  placed ", by.size(), " background apertures");

        // Create aperture mask
        vec2d aper_mask; {
            uint_t y0 = det_img.dims[0]/2;
            uint_t x0 = det_img.dims[1]/2;
            aper_mask = circular_mask(det_img.dims, bg_radius, y0, x0);

            // Write region
            write_ds9_region(outdir+img.base_filename+"_aper.reg",
                det_wcs, vec1u({y0}), vec1u({x0}), bg_radius);
        }

        if (verbose) {
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
        rdata[where(img.mask == 0.0)] = 0;

        img.flux_bg = replicate(dnan, bx.size());
        vec1d bg_area(bx.size());
        for (uint_t i : range(bx)) {
            // Make background aperture mask
            vec2d bg_mask = circular_mask(aper_mask.dims, bg_radius, by[i], bx[i]);
            bg_area[i] = total(bg_mask);
            // Compute fractional area covered inside the mask
            double covcor = bg_area[i]/total(bg_mask*img.mask);
            // Sum pixels and correct for coverage
            img.flux_bg[i] = total(rdata*bg_mask*img.mask)*covcor;
        }

        // Clip strong background outliers
        vec1u idgbg = where(sigma_clip(img.flux_bg, 5.0));
        if (idgbg.size() < 3) {
            // Not enough points left, be less stringent
            idgbg = where(sigma_clip(img.flux_bg, 3.0));
            if (idgbg.size() < 3) {
                // Still not enough, just use everything
                idgbg = uindgen(img.flux_bg.size());
            }
        }

        // Compute and subtract background level
        img.background = median(img.flux_bg[idgbg]/bg_area[idgbg]);
        rdata -= img.background;

        // Compute aperture flux
        img.flux = total(rdata*aper_mask*img.mask);

        // Compute uncertainty from background apertures
        // (the second term accounts for the uncertainty on the background level)
        img.flux_err = stddev(img.flux_bg[idgbg])*sqrt(1.0 + 1.0/idgbg.size());

        // Compute aperture correction from this image's PSF (assumes point source!)
        vec2d model = translate(rpsf, y0 - rpsf.dims[0]/2, x0 - rpsf.dims[1]/2);
        img.apcor = 1.0/total(model*aper_mask*img.mask);

        // Apply aperture correction
        img.flux *= img.apcor;
        img.flux_err *= img.apcor;

        // Write regridded, convolved and background subtracted image
        rdata[where(img.mask == 0.0)] = dnan;
        fits::write(outdir+img.base_filename+"_regrid.fits",
            rdata, img.regridded_hdr
        );
    }

    for (auto& img : images) {
        if (verbose) note("flux for ", img.filename, ": ", img.flux, " +/- ", img.flux_err);
    }

    if (verbose) {
        note("saving catalog...");
    }

    vec1s imgfile;
    vec1f flux, flux_err, apcor, background;
    vec1s bands; vec1u eazy_bands; vec1f lambda;

    for (auto& img : images) {
        imgfile.push_back(img.filename);
        flux.push_back(img.flux);
        flux_err.push_back(img.flux_err);
        background.push_back(img.background);
        apcor.push_back(img.apcor);
        bands.push_back(img.band);
        lambda.push_back(img.lambda);
        eazy_bands.push_back(img.eazy_band);
    }

    fits::write_table(outdir+"fluxes.fits", ftable(
        imgfile, apcor, background, flux, flux_err, bands, lambda, eazy_bands
    ));

    return 0;
}
