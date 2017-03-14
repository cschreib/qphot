#include <phypp.hpp>
#include "image.hpp"
#include "ds9.hpp"

int phypp_main(int argc, char* argv[]) {
    std::string outdir = "./";
    bool verbose = false;
    std::string exclude;
    std::string include;
    std::string fit_mask;
    std::string method = "radial";
    std::string template_file;
    double ra = dnan, dec = dnan;
    double axis_ratio = 1.0;
    double angle = 0.0;

    read_args(argc-1, argv+1, arg_list(
        outdir, verbose, exclude, include, ra, dec, fit_mask, method, axis_ratio, angle,
        template_file
    ));

    outdir = file::directorize(outdir);
    file::mkdir(outdir);

    if (!is_finite(ra) || !is_finite(dec)) {
        error("please provide ra=... and dec=... of the star to subtract");
        return 1;
    }

    // Get files and their properties
    vec<1,image_source_t> images;
    if (!read_image_source_list(argv[1], images)) {
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

    vec1d tpl_r, tpl_f;
    if (!template_file.empty()) {
        fits::read_table(template_file, "d", tpl_r, "f", tpl_f);
    }

    for (auto& img : images) {
        fits::input_image iimg(img.filename);

        if (!iimg.is_image()) {
            error("'", img.filename, "' is not a FITS image");
            return 1;
        }

        fits::header hdr = iimg.read_header();
        astro::wcs w = astro::wcs(hdr);
        double aspix;
        if (!astro::get_pixel_size(w, aspix)) {
            error("could not find pixel size in '", img.filename, "'");
            return 1;
        }

        vec2d data;
        iimg.read(data);

        vec2d psf;
        if (!img.psffile.empty()) {
            // Read the PSF
            fits::read(img.psffile, psf);

            // Adapt the PSF image to the dimensions of the analyzed image
            // and make sure it is centered
            vec1i idm = mult_ids(psf, max_id(psf));
            psf = recenter(psf, idm[0], idm[1], psf.dims);
            psf[where(!is_finite(psf))] = 0;
        } else {
            // Create a PSF if it is not provided
            psf = gaussian_profile(data.dims, img.seeing/2.355/aspix,
                data.dims[0]/2, data.dims[1]/2
            );
        }

        std::string bad_file = file::remove_extension(img.filename)+"-bad.reg";
        if (!file::exists(bad_file)) {
            bad_file = file::remove_extension(img.filename)+"_bad.reg";
            if (!file::exists(bad_file)) {
                bad_file = "";
            }
        }

        if (!bad_file.empty()) {
            // Mask bad pixels
            vec2d bad;
            read_ds9_region_circles_physical(bad_file, w, bad);
            vec2d bad_mask(data.dims);
            for (uint_t i : range(bad.dims[0])) {
                bad_mask += circular_mask(data.dims, bad(i,2), bad(i,1), bad(i,0));
            }

            data[where(bad_mask > 0.5)] = dnan;
        }

        double fzero = fraction_of(data == 0);
        double fnfin = fraction_of(!is_finite(data));
        if (fzero > 0.4) {
            warning("ignoring image '", img.filename, "' because it contains ",
                round(fzero*100), "% of zero values");
            continue;
        } else if (fnfin > 0.4) {
            warning("ignoring image '", img.filename, "' because it contains ",
                round(fnfin*100), "% of invalid values");
            continue;
        }

        vec2b mask = !is_finite(data);
        vec2d maskreg;
        if (file::exists(fit_mask)) {
            read_ds9_region_circles_physical(fit_mask, w, maskreg);
            vec2d bad_mask(data.dims);
            for (uint_t i : range(maskreg.dims[0])) {
                bad_mask += circular_mask(data.dims, maskreg(i,2), maskreg(i,1), maskreg(i,0));
            }

            mask = mask || bad_mask > 0.5;
        }

        double sx, sy;
        astro::ad2xy(w, ra, dec, sx, sy);
        sx -= 1.0; sy -= 1.0;

        if (method == "radial") {
            vec2d d = generate_img(data.dims, [&](double y, double x) {
                double ody = y - sy;
                double odx = x - sx;
                double dx = cos(angle*dpi/180)*odx + sin(angle*dpi/180)*ody;
                double dy = cos(angle*dpi/180)*ody - sin(angle*dpi/180)*odx;
                dx *= axis_ratio;
                return sqrt(sqr(dy) + sqr(dx));
            });

            vec2d db = make_bins(0.0, max(d), uint_t(ceil(max(d)/3.00)));
            vec1d dx = bin_center(db);

            vec2d tdata = data;
            tdata[where(mask)] = dnan;
            vec1d df(dx.dims);

            histogram(d, db, [&](uint_t i, vec1u ids) {
                if (ids.empty()) {
                    df[i] = dnan;
                } else {
                    df[i] = median(tdata[ids]);
                }
            });

            vec1u idg = where(is_finite(df));
            if (!idg.empty()) {
                data -= reform(interpolate(df[idg], dx[idg], flatten(d)), data.dims);
                fits::write_table(file::remove_extension(img.filename)+"_profile.fits",
                    "d", dx[idg]*aspix, "f", df[idg]
                );
            } else {
                warning("could not clean ", img.filename);
            }
        } else if (method == "psf") {
            linfit_result bres;
            vec2d bpsf;
            bres.chi2 = dinf;
            int_t dx_max = 5, dy_max = 5;
            double dp = 0.2;
            int_t np = psf.dims[0]/2;
            int_t ni = data.dims[0]/2;
            vec1u idf = where(!mask);

            sx -= ni;
            sy -= ni;

            for (int_t dy = -dy_max; dy <= dy_max; ++dy)
            for (int_t dx = -dx_max; dx <= dx_max; ++dx) {
                vec2d npsf = translate(psf, sy+dy*dp, sx+dx*dp);
                npsf = npsf((np-ni)-_-(np+ni),(np-ni)-_-(np+ni));

                auto res = linfit(data[idf], 1.0, 1.0, npsf[idf]);
                if (res.success && res.chi2 < bres.chi2) {
                    bres = res;
                    bpsf = npsf;
                }
            }

            if (bres.success) {
                data -= bpsf*bres.params[1];
            } else {
                warning("could not clean ", img.filename);
            }
        } else if (method == "template") {
            vec2d d = generate_img(data.dims, [&](double y, double x) {
                double ody = y - sy;
                double odx = x - sx;
                double dx = cos(angle*dpi/180)*odx + sin(angle*dpi/180)*ody;
                double dy = cos(angle*dpi/180)*ody - sin(angle*dpi/180)*odx;
                dx *= axis_ratio;
                return aspix*sqrt(sqr(dy) + sqr(dx));
            });

            vec2d model = reform(interpolate(tpl_f, tpl_r, flatten(d)), data.dims);

            linfit_result bres;
            vec2d bpsf;
            bres.chi2 = dinf;
            int_t dx_max = 5, dy_max = 5;
            double dp = 0.2;
            vec1u idf = where(!mask);

            for (int_t dy = -dy_max; dy <= dy_max; ++dy)
            for (int_t dx = -dx_max; dx <= dx_max; ++dx) {
                vec2d npsf = convolve2d(translate(model, dy*dp, dx*dp), psf);
                auto res = linfit(data[idf], 1.0, 1.0, npsf[idf]);
                if (res.success && res.chi2 < bres.chi2) {
                    bres = res;
                    bpsf = npsf;
                }
            }

            if (bres.success) {
                data -= bpsf*bres.params[1];
            } else {
                warning("could not clean ", img.filename);
            }
        }

        fits::output_image oimg(outdir+file::get_basename(img.filename));
        oimg.write(data);
        oimg.write_header(hdr);
    }

    return 0;
}
