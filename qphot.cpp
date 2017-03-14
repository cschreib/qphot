#include <phypp.hpp>
#include "state.hpp"

int phypp_main(int argc, char* argv[]) {
    options_t opts;

    read_args(argc-1, argv+1, arg_list(
        opts.aperture, opts.outdir, opts.verbose, opts.exclude, opts.include, opts.det_seeing,
        opts.error_file, opts.min_bg_aper, opts.clean_model, opts.reuse, opts.clean_ratio,
        opts.clean_threshold, opts.clean_conv_seeing, opts.clean_group, opts.model_dance,
        opts.dance_step, opts.dance_max, opts.dance_chi2_range, opts.det_min_area,
        opts.debug_clean, opts.no_neighbor_mask_background, opts.save_models, opts.det_threshold,
        opts.bg_threshold, opts.hri_min_snr, opts.extra, opts.nocuts, opts.priors, opts.cutout_size,
        opts.show_progress, opts.save_interval, opts.threads
    ));

    if (opts.clean_model != "psf" && opts.clean_model != "hri" &&
        opts.clean_model != "nearest_hri" && opts.clean_model != "clean") {
        error("unknown clean model '", opts.clean_model, "'");
        note("allowed values: psf, hri, nearest_hri, clean");
        return 1;
    }

    if (opts.clean_group && opts.clean_model == "nearest_hri") {
        warning("cannot use group cleaning with 'nearest_hri' cleaning");
        opts.clean_group = false;
    }

    if (opts.save_models && opts.nocuts) {
        warning("disabling 'save_models' option because 'nocuts' is set");
        opts.save_models = false;
    }

    if (opts.show_progress && opts.verbose) {
        warning("disabling 'show_progress' option because 'verbose' is set");
        opts.show_progress = false;
    }

    opts.outdir = file::directorize(opts.outdir);
    file::mkdir(opts.outdir);

    // Get files and their properties
    if (opts.verbose) {
        note("reading image list...");
    }
    vec<1,image_source_t> image_sources;
    if (!read_image_source_list(argv[1], image_sources)) {
        return 1;
    }

    if (!opts.extra.empty()) {
        if (!read_image_source_list(opts.extra, image_sources)) {
            return 1;
        }
    }

    // Remove bands we don't want
    if (!opts.exclude.empty() || !opts.include.empty()) {
        vec1u idb;
        for (uint_t i : range(image_sources)) {
            if ((!opts.exclude.empty() && regex_match(image_sources[i].filename, opts.exclude)) ||
                (!opts.include.empty() && !regex_match(image_sources[i].filename, opts.include))) {
                idb.push_back(i);
            }
        }

        inplace_remove(image_sources, idb);
    }

    // Check validity and initialize
    if (opts.verbose) {
        note("checking input images...");
    }

    if (!check_image_source_list(image_sources)) {
        return 1;
    }

    if (!opts.priors.empty()) {
        // Read a set of positions from a catalog and extract cutouts on the fly
        vec1s sid;
        vec1d ra, dec;
        fits::read_table(opts.priors, ftable(sid, ra, dec));

        vec2f flux, flux_err, apcor, background;
        vec2u num_bg;

        vec1s imgfile;
        vec1s bands;
        vec1u eazy_bands;
        vec1f lambda;

        for (auto& img : image_sources) {
            imgfile.push_back(img.filename);
            bands.push_back(img.band);
            lambda.push_back(img.lambda);
            eazy_bands.push_back(img.eazy_band);
        }

        auto save_catalog = [&]() {
            if (opts.verbose) {
                note("saving catalog...");
            }

            fits::write_table(opts.outdir+"fluxes.fits", ftable(
                sid, ra, dec, flux, flux_err, apcor, background, num_bg,
                imgfile, bands, lambda, eazy_bands
            ));
        };

        if (!opts.error_file.empty()) {
            file::remove(opts.error_file);
        }

        bool saved = false;
        auto pg = progress_start(sid.size());
        for (uint_t i : range(sid)) {
            state_t st(image_sources, opts);

            st.outdir = opts.outdir+sid[i]+"/";
            if (!opts.nocuts) {
                file::mkdir(st.outdir);
            }

            if (!opts.error_file.empty()) {
                st.errfile.open(opts.error_file, std::ios::app);
            }

            st.extract_cutouts(ra[i], dec[i]);
            st.read_psfs();
            st.build_detection_image();
            st.clean_images();
            st.extract_fluxes();

            vec1f tflx, terr, tbg, tac;
            vec1u tnbg;
            for (auto& img : st.images) {
                tflx.push_back(img.flux);
                terr.push_back(img.flux_err);
                tbg.push_back(img.background);
                tac.push_back(img.apcor);
                tnbg.push_back(img.num_bg);
            }

            append<0>(flux,       reform(std::move(tflx), 1, tflx.dims));
            append<0>(flux_err,   reform(std::move(terr), 1, terr.dims));
            append<0>(apcor,      reform(std::move(tac),  1, tac.dims));
            append<0>(background, reform(std::move(tbg),  1, tbg.dims));
            append<0>(num_bg,     reform(std::move(tnbg), 1, tnbg.dims));
            saved = false;

            if (opts.save_interval != 0 && i % opts.save_interval == 0 && i != 0) {
                save_catalog();
                saved = true;
                return 0;
            }

            if (opts.show_progress) progress(pg);
        }

        if (!saved) {
            save_catalog();
        }
    } else {
        for (auto& img : image_sources) {
            if (img.input_images.size() > 1) {
                error("secfits files not supported in single source mode");
                error("please provide only FITS images");
                return 1;
            }
        }

        // Extract a source at the center of the provided images
        state_t st(image_sources, opts);

        if (!opts.error_file.empty()) {
            st.errfile.open(opts.error_file);
        }

        st.read_cutouts();
        st.read_psfs();
        st.build_detection_image();
        st.clean_images();
        st.extract_fluxes();

        if (opts.verbose) {
            for (auto& img : st.images) {
                note("flux for ", img.source.filename, ": ", img.flux, " +/- ", img.flux_err);
            }
        }

        if (opts.verbose) {
            note("saving catalog...");
        }

        vec1s imgfile;
        vec1f flux, flux_err, apcor, background;
        vec1s bands; vec1u eazy_bands; vec1f lambda;
        vec1u num_bg;

        for (auto& img : st.images) {
            flux.push_back(img.flux);
            flux_err.push_back(img.flux_err);
            background.push_back(img.background);
            apcor.push_back(img.apcor);
            num_bg.push_back(img.num_bg);
            imgfile.push_back(img.source.filename);
            bands.push_back(img.source.band);
            lambda.push_back(img.source.lambda);
            eazy_bands.push_back(img.source.eazy_band);
        }

        fits::write_table(opts.outdir+"fluxes.fits", ftable(
            imgfile, apcor, background, flux, flux_err, bands, lambda, eazy_bands, num_bg
        ));
    }

    return 0;
}
