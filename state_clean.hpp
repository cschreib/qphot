void state_t::clean_images() {
    if (!good || !has_clean) return;

    if (opts.verbose) {
        note("cleaning images...");
    }

    // Clean homogenized images
    if (has_homogenize && opts.clean_group) {
        if (opts.verbose) note("cleaning homogenized images...");

        // Make sure the images are homogenized first
        vec2b fitmask = replicate(true, det_img.dims);
        uint_t nhomo = 0;
        for (auto& img : images) {
            if (!img.covered || !img.source.clean_neighbors || img.source.noconvolve) continue;

            if (opts.verbose) note("  preparing ", img.source.short_name);

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

            fitmask = fitmask && is_finite(rdata);
            img.mask = vec2d(is_finite(rdata));
            rdata[where(img.mask < 0.5)] = 0;
            ++nhomo;
        }

        // Build models
        if (opts.verbose) note("  building models...");
        vec2d models(1 + segments.id.size(), det_img.dims[0]*det_img.dims[1]);
        models(0,_) = 1.0; // background
        vec1u gmodel; gmodel.push_back(0);

        if (opts.clean_model == "psf") {
            // Assume point sources
            for (uint_t s : range(segments.id)) {
                models(s+1,_) = flatten(make_point_source_integer(
                    det_psf, det_img.dims, segments.py[s], segments.px[s]
                ));
            }
        } else if (opts.clean_model == "clean") {
            uint_t imod = 1;
            auto c2d = batch_convolve2d(det_psf);

            foreach_segment(seg, segments.origin, [&](uint_t is, const vec1u& ids) {
                // Get pixels of that segment
                vec2d tmod(det_img.dims);
                tmod[ids] = det_model[ids]/max(det_model[ids]);
                // Convolve to the resolution of the fitted image
                // NB: we have to assume the same PSF for all homogenized images
                tmod = c2d.convolve(tmod);
                // Save model
                models(is+imod,_) = flatten(tmod);
            });

            imod += segments.id.size();
        } else if (opts.clean_model == "hri") {
            auto& hdata = hri[0]->regridded_data_noconvol;
            if (hdata.empty()) {
                // Regrid, if not done at the detection stage already
                if (opts.verbose) note("  regridding HRI...");
                regrid(*hri[0]);
            }

            double rms = 1.48*mad(hdata);
            uint_t imod = 1;
            auto c2d = batch_convolve2d(det_psf);

            foreach_segment(seg, segments.origin, [&](uint_t is, const vec1u& ids) {
                // Get pixels of that segment
                vec2d tmod(det_img.dims);
                tmod[ids] = hdata[ids];
                double tt = total(tmod[ids])/(rms*sqrt(ids.size()));
                if (tt > opts.hri_min_snr) {
                    // Keep the HRI model
                    tmod /= max(hdata[ids]);
                    // Convolve to the resolution of the fitted image
                    tmod = c2d.convolve(tmod);
                    models(is+imod,_) = flatten(tmod);
                } else {
                    // Not enough SNR, make that a point source
                    models(is+imod,_) = flatten(make_point_source_integer(
                        det_psf, det_img.dims, segments.py[is], segments.px[is]
                    ));
                }
            });

            imod += segments.id.size();
        } else if (opts.clean_model == "nearest_hri") {
            write_error("trying to clean with 'nearest_hri' with group cleaning enabled");
            good = false;
            return;
        }

        // Inspect models and flag out the ones that are zero
        for (uint_t s : range(segments.id)) {
            if (total(sqr(models(s+1,_))) > 0) {
                gmodel.push_back(s+1);
            } else {
                if (s == isource) {
                    write_error("fit model for the target only contains zero values");
                    good = false;
                    return;
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
            good = false;
            return;
        }

        int_t dance_step_pix = max(1, ceil(opts.dance_step/det_aspix));
        int_t dance_nstep = ceil(opts.dance_max/det_aspix/dance_step_pix);
        auto dofit = [&](int_t dy, int_t dx) {
            uint_t ihomo = 0;
            for (auto& img : images) {
                if (!img.covered || !img.source.clean_neighbors || img.source.noconvolve) continue;
                if (img.source.nodance && (dy != 0 || dx != 0)) continue;

                fit_neighbors(img, batch, models, idf, gmodel, dance_step_pix*dy, dance_step_pix*dx);

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

        if (opts.verbose) note("  fitting...");

        // Do the fit
        dofit(0, 0);

        if (opts.model_dance && dance_nstep > 0) {
            if (opts.verbose) {
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
                if (!img.covered || !img.source.clean_neighbors || img.source.noconvolve) continue;

                if (fit_dx[ihomo] != 0 || fit_dy[ihomo] != 0) {
                    // Shift the image for flux measurement later on
                    img.regridded_data = translate_integer(img.regridded_data,
                        -dance_step_pix*fit_dy[ihomo], -dance_step_pix*fit_dx[ihomo]
                    );
                    // Don't forget to shift the mask too!
                    img.mask = translate_integer(img.mask,
                        -dance_step_pix*fit_dy[ihomo], -dance_step_pix*fit_dx[ihomo]
                    );
                }

                fits::setkey(img.regridded_hdr, "DANCEX",
                    img.regridded_aspix*dance_step_pix*fit_dx[ihomo], "dance offset [arcsec]");
                fits::setkey(img.regridded_hdr, "DANCEY",
                    img.regridded_aspix*dance_step_pix*fit_dy[ihomo], "dance offset [arcsec]");
                fits::setkey(img.regridded_hdr, "DANCES",
                    img.regridded_aspix*dance_step_pix, "dance step [arcsec]");
                fits::setkey(img.regridded_hdr, "DANCER",
                    img.regridded_aspix*dance_nstep*dance_step_pix, "dance radius [arcsec]");

                ++ihomo;
            }
        }

        if (opts.verbose) note("  subtracting neighbors...");

        // Remove all the fitted segments from the image
        uint_t ihomo = 0;
        for (auto& img : images) {
            if (!img.covered || !img.source.clean_neighbors || img.source.noconvolve) continue;

            auto& rdata = img.regridded_data;

            for (uint_t m : range(gmodel)) {
                if (gmodel[m] == isource+1) {
                    // Keep our target there
                    continue;
                }

                if (!is_finite(fit_fluxes(ihomo,m))) {
                    write_warning("in the fit for image ", img.source.short_name);
                    if (gmodel[m] == 0) {
                        write_warning("the fitted value of the background was invalid");
                        write_warning("something is very wrong!");
                    } else {
                        write_warning("the fitted value of source ", gmodel[m]-1, " was invalid");
                    }
                    continue;
                }

                vec2d tmod = reform(models(gmodel[m],_), rdata.dims);
                rdata -= fit_fluxes(ihomo,m)*tmod;
            }

            rdata[where(img.mask < 0.5)] = dnan;
            img.cleaned = true;
            ++ihomo;
        }
    }

    for (auto& img : images) {
        if (!img.covered || !img.source.clean_neighbors || img.cleaned) continue;

        if (opts.verbose) note("cleaning ", img.source.short_name);

        auto& rdata = img.regridded_data;
        auto& rpsf = img.regridded_psf;

        if (rdata.empty()) {
            // Regrid, if not done at the detection stage already
            if (opts.verbose) note("  regridding...");
            regrid(img);
        }

        if (!img.source.noconvolve && !img.convolved) {
            if (opts.verbose) note("  convolving for homogenization...");
            convolve_to(img, worst_seeing);
        }

        if (opts.verbose) note("  building models...");

        // Build models
        vec2d models(1 + segments.id.size(), rdata.dims[0]*rdata.dims[1]);
        models(0,_) = 1.0; // background
        vec1u gmodel; gmodel.push_back(0);

        if (opts.clean_model == "psf") {
            // Assume point sources
            for (uint_t s : range(segments.id)) {
                models(s+1,_) = flatten(make_point_source_integer(
                    rpsf, rdata.dims, segments.py[s], segments.px[s]
                ));
            }
        } else if (opts.clean_model == "clean") {
            uint_t imod = 1;
            auto c2d = batch_convolve2d(rpsf);

            foreach_segment(seg, segments.origin, [&](uint_t is, const vec1u& ids) {
                // Get pixels of that segment
                vec2d tmod(det_img.dims);
                tmod[ids] = det_model[ids]/max(det_model[ids]);
                // Convolve to the resolution of the fitted image
                tmod = c2d.convolve(tmod);
                // Save model
                models(is+imod,_) = flatten(tmod);
            });

            imod += segments.id.size();
        } else if (opts.clean_model == "hri") {
            auto& hdata = hri[0]->regridded_data_noconvol;
            if (hdata.empty()) {
                // Regrid, if not done at the detection stage already
                if (opts.verbose) note("  regridding HRI...");
                regrid(*hri[0]);
            }

            double rms = 1.48*mad(hdata);
            uint_t imod = 1;
            auto c2d = batch_convolve2d(rpsf);

            foreach_segment(seg, segments.origin, [&](uint_t is, const vec1u& ids) {
                // Get pixels of that segment
                vec2d tmod(det_img.dims);
                tmod[ids] = hdata[ids];
                double tt = total(tmod[ids])/(rms*sqrt(ids.size()));
                if (tt > opts.hri_min_snr) {
                    // Keep the HRI model
                    tmod /= max(tmod[ids]);
                    // Convolve to the resolution of the fitted image
                    tmod = c2d.convolve(tmod);
                    models(is+imod,_) = flatten(tmod);
                } else {
                    // Not enough SNR, make that a point source
                    models(is+imod,_) = flatten(make_point_source_integer(
                        rpsf, rdata.dims, segments.py[is], segments.px[is]
                    ));
                }
            });

            imod += segments.id.size();
        } else if (opts.clean_model == "nearest_hri") {
            vec1d dlam(hri.size());
            vec1d hrms(hri.size());
            for (uint_t h : range(hri)) {
                if (hri[h]->regridded_data_noconvol.empty()) {
                    // Regrid, if not done at the detection stage already
                    if (opts.verbose) note("  regridding HRI...");
                    regrid(*hri[h]);
                }

                dlam[h] = log10(img.source.lambda/hri[h]->source.lambda);
                hrms[h] = 1.48*mad(hri[h]->regridded_data_noconvol);
            }

            vec1u shri = sort(dlam);
            auto c2d = batch_convolve2d(rpsf);

            uint_t imod = 1;
            foreach_segment(seg, segments.origin, [&](uint_t is, const vec1u& ids) {
                // Get pixels of that segment
                vec2d tmod(det_img.dims);

                uint_t ihri = 0;
                do {
                    tmod[ids] = hri[shri[ihri]]->regridded_data_noconvol[ids];
                    double tt = total(tmod[ids])/(hrms[shri[ihri]]*sqrt(ids.size()));
                    if (tt > opts.hri_min_snr) {
                        // Keep the HRI model
                        tmod /= max(tmod[ids]);
                        // Convolve to the resolution of the fitted image
                        tmod = c2d.convolve(tmod);
                        models(is+imod,_) = flatten(tmod);
                        break;
                    } else {
                        // Not enough SNR, try next HRI
                        ++ihri;
                    }
                } while (ihri < hri.size());

                if (ihri == hri.size()) {
                    // Not enough SNR, make that a point source
                    models(is+imod,_) = flatten(make_point_source_integer(
                        rpsf, rdata.dims, segments.py[is], segments.px[is]
                    ));
                }
            });

            imod += segments.id.size();
        }

        // Inspect models and flag out the ones that are zero
        for (uint_t s : range(segments.id)) {
            if (total(sqr(models(s+1,_))) > 0) {
                gmodel.push_back(s+1);
            } else {
                if (s == isource) {
                    write_error("fit model for the target only contains zero values");
                    good = false;
                    return;
                } else {
                    write_warning("model for source ", s, " is zero");
                    write_warning("this source will be ignored in the fit");
                }
            }
        }

        if (opts.verbose) note("  fitting...");

        // Do the fit
        img.mask = vec2d(is_finite(rdata));
        rdata[where(img.mask < 0.5)] = 0;
        vec1u idf = img.idf;
        vec1d err = replicate(1.0, idf.size());
        vec1d fit_fluxes(gmodel.size());
        int_t fit_dx = 0, fit_dy = 0;
        float fit_chi2 = finf;
        auto batch = linfit_pack_batch<vec1d>(err, models(gmodel,idf));
        if (!batch.fr.success) {
            write_error("could not perform the fit to remove the neighbors");
            good = false;
            return;
        }

        int_t dance_step_pix = max(1, ceil(opts.dance_step/img.regridded_aspix));
        int_t dance_nstep = ceil(opts.dance_max/img.regridded_aspix/dance_step_pix);
        auto dofit = [&](int_t dy, int_t dx) {
            fit_neighbors(img, batch, models, idf, gmodel, dance_step_pix*dy, dance_step_pix*dx);

            if (batch.fr.success && batch.fr.chi2 < fit_chi2) {
                fit_chi2 = batch.fr.chi2;
                fit_fluxes = batch.fr.params;
                fit_dy = dy;
                fit_dx = dx;
            }
        };

        dofit(0, 0);

        // Adjust positions if asked
        if (opts.model_dance && dance_nstep > 0 && !img.source.nodance) {
            if (opts.verbose) {
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
                // Don't forget to shift the mask too!
                img.mask = translate_integer(img.mask,
                    -dance_step_pix*fit_dy, -dance_step_pix*fit_dx
                );
            }

            fits::setkey(img.regridded_hdr, "DANCEX",
                img.regridded_aspix*dance_step_pix*fit_dx, " [arcsec]");
            fits::setkey(img.regridded_hdr, "DANCEY",
                img.regridded_aspix*dance_step_pix*fit_dy, " [arcsec]");
            fits::setkey(img.regridded_hdr, "DANCES",
                img.regridded_aspix*dance_step_pix, " [arcsec]");
            fits::setkey(img.regridded_hdr, "DANCER",
                img.regridded_aspix*dance_nstep*dance_step_pix, " [arcsec]");
        }

        if (opts.verbose) note("  subtracting neighbors...");

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
                good = false;
                return;
            }

            vec2d tmod = reform(models(gmodel[m],_), rdata.dims);
            rdata -= fit_fluxes[m]*tmod;
        }

        if (opts.save_models) {
            vec2d tmodel(rdata.dims);
            for (uint_t m : range(gmodel)) {
                if (!is_finite(fit_fluxes[m])) continue;

                vec2d tmod = reform(models(gmodel[m],_), rdata.dims);
                tmodel += fit_fluxes[m]*tmod;
            }

            fits::output_image oimg(outdir+img.source.short_name+"_regrid_model.fits");
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
                id[im] = segments.id[gmodel[m]-1];
                sx[im] = segments.px[gmodel[m]-1];
                sy[im] = segments.py[gmodel[m]-1];

                ++im;
            }

            vec1u idsource = gmodel[where(gmodel != 0)];

            fits::output_table otbl(outdir+img.source.short_name+"_regrid_model_catalog.fits");
            otbl.reach_hdu(1);
            otbl.write_columns("id", id, "x", sx, "y", sy, "flux", flux);
        }

        img.regridded_data[where(img.mask < 0.5)] = dnan;
        img.cleaned = true;
    }
}
