#include <phypp.hpp>
#include "state.hpp"

int phypp_main(int argc, char* argv[]) {
    astro::wcs::enable_errors();

    options_t opts;

    read_args(argc-1, argv+1, arg_list(
        opts.aperture, opts.outdir, opts.verbose, opts.exclude, opts.include, opts.det_seeing,
        opts.error_file, opts.min_bg_aper, opts.clean_model, opts.reuse, opts.reuse_cutouts,
        opts.clean_ratio, opts.clean_threshold, opts.clean_conv_seeing, opts.clean_group,
        opts.model_dance, opts.dance_step, opts.dance_max, opts.dance_chi2_range, opts.det_min_area,
        opts.debug_clean, opts.no_neighbor_mask_background, opts.save_models, opts.det_threshold,
        opts.bg_threshold, opts.hri_min_snr, opts.extra, opts.nocuts, opts.priors, opts.cutout_size,
        opts.show_progress, opts.save_interval, opts.threads, opts.merge, opts.update,
        opts.search_radius, opts.prior_ra, opts.prior_dec, opts.only_sources
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

    if (opts.threads > 1 && fits_is_reentrant() == 0) {
        warning("the CFITSIO library does not support multithreading, will run in single thread");
        opts.threads = 1;
    }
    if (opts.threads == 0) {
        opts.threads = 1;
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


    vec1s osid;
    vec1d ora, odec;
    vec1u oid;

    if (!opts.priors.empty()) {
        fits::read_table(opts.priors, "sid", osid, "ra", ora, "dec", odec);

        oid = replicate(npos, osid.size());

        if (!opts.only_sources.empty()) {
            vec1u id1, id2;
            match(osid, opts.only_sources, id1, id2);

            osid = osid[id1];
            ora = ora[id1];
            odec = odec[id1];
            oid = oid[id1];
        }
    } else if (is_finite(opts.prior_ra) && is_finite(opts.prior_dec)) {
        osid = {""};
        ora = {opts.prior_ra};
        odec = {opts.prior_dec};
        oid = {npos};
    }

    if (!osid.empty()) {
        // Extract cutouts on the fly
        vec1s sid;
        vec1d ra, dec;
        vec2f flux, flux_err, apcor, background;
        vec2u num_bg;

        vec1s imgfile;
        vec1s bands;
        vec1u eazy_bands;
        vec1f lambda;

        if ((opts.merge || opts.update) && file::exists(opts.outdir+"fluxes.fits")) {
            fits::read_table(opts.outdir+"fluxes.fits", ftable(
                sid, ra, dec, flux, flux_err, apcor, background, num_bg,
                imgfile, bands, lambda, eazy_bands
            ));

            if (opts.verbose) {
                note("found ", sid.size(), " sources in the existing catalog");
            }

            vec1u id1, id2;
            match(osid, sid, id1, id2);

            if (opts.merge) {
                inplace_remove(osid, id1);
                inplace_remove(ora, id1);
                inplace_remove(odec, id1);
                inplace_remove(oid, id1);

                if (opts.verbose) {
                    note("will analyze ", osid.size(), " new sources");
                }
            } else if (opts.update) {
                oid[id1] = id2;

                if (opts.verbose) {
                    if (osid.size() - id1.size() > 0 && !id1.empty()) {
                        note("will analyze ", osid.size() - id1.size(), " new sources and update ",
                            id1.size(), " sources");
                    } else if (id1.empty()) {
                        note("will analyze ", osid.size(), " new sources");
                    } else {
                        note("will update ", id1.size(), " sources");
                    }
                }
            }
        } else {
            if (opts.verbose) {
                note("will analyze ", osid.size(), " sources");
            }

            for (auto& img : image_sources) {
                imgfile.push_back(img.filename);
                bands.push_back(img.band);
                lambda.push_back(img.lambda);
                eazy_bands.push_back(img.eazy_band);
            }
        }

        if (ora.empty()) {
            if (opts.verbose) {
                note("no source to analyze");
            }

            return 0;
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

        struct output_t {
            uint_t id;
            vec1d flux, flux_err, background, apcor;
            vec1u num_bg;
        };

        // Thread safe processing function
        auto process_source = [&](uint_t i) {
            output_t o;
            o.id = i;

            state_t st(image_sources, opts);

            if (osid[i] != "") {
                st.outdir = opts.outdir+osid[i]+"/";
            } else {
                st.outdir = opts.outdir;
            }

            if (!opts.nocuts) {
                file::mkdir(st.outdir);
            }

            if (!opts.error_file.empty()) {
                st.errfile.open(opts.error_file, std::ios::app);
            }

            st.extract_cutouts(ora[i], odec[i]);
            st.read_psfs();
            st.build_detection_image();
            st.apply_filters();
            st.clean_images();
            st.extract_fluxes();

            for (auto& img : st.images) {
                o.flux.push_back(img.flux);
                o.flux_err.push_back(img.flux_err);
                o.background.push_back(img.background);
                o.apcor.push_back(img.apcor);
                o.num_bg.push_back(img.num_bg);
            }

            return o;
        };

        // Non-thread safe function to compile results
        bool saved = true;
        auto save_source = [&](output_t o) {
            uint_t nimg = o.flux.size();

            if (oid[o.id] == npos) {
                // New source
                sid.push_back(osid[o.id]);
                ra.push_back(ora[o.id]);
                dec.push_back(odec[o.id]);
                append<0>(flux,       reform(std::move(o.flux),       1,     nimg));
                append<0>(flux_err,   reform(std::move(o.flux_err),   1,     nimg));
                append<0>(apcor,      reform(std::move(o.apcor),      1,     nimg));
                append<0>(background, reform(std::move(o.background), 1,     nimg));
                append<0>(num_bg,     reform(std::move(o.num_bg),     1,     nimg));
            } else {
                // Update source
                uint_t j = oid[o.id];
                sid[j] = osid[o.id];
                ra[j] = ora[o.id];
                dec[j] = odec[o.id];
                flux(j,_) = o.flux;
                flux_err(j,_) = o.flux_err;
                apcor(j,_) = o.apcor;
                background(j,_) = o.background;
                num_bg(j,_) = o.num_bg;
            }

            saved = false;

            if (opts.save_interval != 0 &&
                sid.size() % opts.save_interval == 0 && sid.size() != 0) {
                save_catalog();
                saved = true;
            }
        };

        if (ora.size() < opts.threads) {
            if (opts.verbose) {
                note("only ", ora.size(), " sources will be extracted, reducing the number of "
                    "threads");
            }

            opts.threads = ora.size();
        }

        if (opts.threads == 1) {
            auto pg = progress_start(osid.size());
            for (uint_t i : range(osid)) {
                if (opts.verbose) note("processing source ", osid[i]);
                save_source(process_source(i));
                if (opts.show_progress) progress(pg);
            }

            if (!saved) {
                save_catalog();
            }
        } else {
            if (opts.threads > 1 && opts.verbose) {
                warning("will run in multiple threads, 'verbose' option disabled");
                opts.verbose = false;
            }

            auto tpool = thread::pool(opts.threads);

            thread::lock_free_queue<output_t> save_queue;
            std::atomic<uint_t> iter(0);

            uint_t ifirst = 0;
            uint_t source_per_thread = osid.size()/opts.threads + 1;
            for (uint_t i : range(opts.threads)) {
                uint_t ilast = min(ifirst + source_per_thread, osid.size());
                vec1u do_ids = uindgen(ilast - ifirst) + ifirst;
                ifirst = ilast;

                tpool[i].start([&,do_ids]() {
                    for (uint_t k : do_ids) {
                        save_queue.push(process_source(k));
                        ++iter;
                    }
                });
            }

            auto pg = progress_start(osid.size());
            if (opts.show_progress) print_progress(pg, iter);

            output_t o;
            while (iter != osid.size()) {
                while (save_queue.pop(o)) {
                    save_source(o);
                    if (opts.show_progress) print_progress(pg, iter);
                }

                thread::sleep_for(0.01);
            }

            for (auto& t : tpool) {
                t.join();
            }

            while (save_queue.pop(o)) {
                save_source(o);
            }

            if (!saved) {
                save_catalog();
            }
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
        st.apply_filters();
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
