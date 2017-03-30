vec2d make_point_source_integer(const vec2d& psf, const std::array<uint_t,2>& dims, int_t y, int_t x) {
    vec2d model(dims);
    int_t y0 = max(0, y - int_t(psf.dims[0]/2));
    int_t y1 = min(int_t(dims[0]), y + int_t(psf.dims[0]/2+1));
    int_t x0 = max(0, x - int_t(psf.dims[1]/2));
    int_t x1 = min(int_t(dims[0]), x + int_t(psf.dims[1]/2+1));

    for (int_t iy = y0; iy < y1; ++iy)
    for (int_t ix = x0; ix < x1; ++ix) {
        int_t ipy = iy - y + int_t(psf.dims[0]/2);
        int_t ipx = ix - x + int_t(psf.dims[1]/2);
        model.safe(uint_t(iy),uint_t(ix)) = psf.safe(uint_t(ipy),uint_t(ipx));
    }

    return model;
}

vec2d make_point_source(const vec2d& psf, const std::array<uint_t,2>& dims, double y, double x) {
    int_t iy = round(y), ix = round(x);
    vec2d npsf = translate(psf, y - iy, x - ix);
    return make_point_source_integer(npsf, dims, iy, ix);
}


void state_t::clean_images() {
    if (!good) return;

    if (has_clean) {
        if (opts.verbose) {
            note("cleaning images...");
        }

        if (opts.verbose) {
            note("cleaning detection image...");
        }

        vec2d model(det_img.dims);
        vec2d det_img_clean; {
            vec2d tdimg = det_img;
            uint_t mid = max_id(det_img);
            uint_t iter = 0;

            vec2d clean_psf = det_psf;
            clean_psf /= max(clean_psf);

            while (tdimg[mid] > opts.clean_threshold) {
                vec1i ids = mult_ids(tdimg, mid);
                vec2d tpsf = make_point_source_integer(clean_psf, tdimg.dims, ids[0], ids[1]);

                if (opts.debug_clean && iter % 20 == 0) {
                    print("clean ", iter, ": max ", tdimg[mid], " pos: ", ids);
                    if (iter % 100 == 0) {
                        fits::write("debug_clean.fits", tdimg);
                    }
                }

                double flx = tdimg[mid]*opts.clean_ratio;
                tdimg -= tpsf*flx;
                model[mid] += flx;

                mid = max_id(tdimg);

                ++iter;
            }

            model[where(model > 0)] += opts.clean_threshold;

            double clean_conv_seeing = opts.clean_conv_seeing;
            if (!is_finite(clean_conv_seeing)) {
                clean_conv_seeing = 0.66*det_seeing;
            }

            vec2d model_psf = gaussian_profile(det_img.dims,
                clean_conv_seeing/2.355/det_aspix, det_img.dims[0]/2, det_img.dims[1]/2
            );

            model_psf /= max(model_psf);

            det_img_clean = convolve2d(model, model_psf);

            model[where(model > 0)] -= opts.clean_threshold;

            if (!opts.nocuts) {
                fits::output_image oimg(outdir+"detection_image_cleaned.fits");
                oimg.reach_hdu(1);
                oimg.write(det_img_clean);
                oimg.write_header(det_hdr);
                oimg.write_keyword("CRATIO", opts.clean_ratio);
                oimg.write_keyword("CTHRESH", opts.clean_threshold);
                oimg.write_keyword("CSEING", clean_conv_seeing);

                oimg.reach_hdu(2);
                oimg.write(model);
                oimg.write_header(det_hdr);
            }
        }

        if (opts.verbose) {
            note("building segmentation map...");
        }

        segment_deblend_params sdp;
        sdp.detect_threshold = opts.det_threshold;
        sdp.deblend_threshold = 0.5*det_seeing/det_aspix;
        sdp.min_area = opts.det_min_area;
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

            if (bestd > sqr(search_radius/det_aspix)) {
                write_warning("the target source was not found in the detection image");
                write_warning("(the closest source was found at ", sqrt(bestd)*det_aspix, " arcsec)");
                write_warning("it will be added as a central point source for the fitting stage");

                isource = sdo.bx.size();
                sdo.id.push_back(max(sdo.id)+1);
                sdo.px.push_back(det_img.dims[1]/2);
                sdo.py.push_back(det_img.dims[0]/2);
                sdo.bx.push_back(sdo.px.back());
                sdo.by.push_back(sdo.py.back());
                sdo.origin.push_back(flat_id(seg, sdo.py.back(), sdo.px.back()));

                if (opts.clean_model == "psf" || opts.clean_model == "clean") {
                    // A point source will suffice
                    sdo.area.push_back(1);
                    sdo.flux.push_back(det_img_clean(sdo.py.back(), sdo.px.back()));
                    model(sdo.py.back(), sdo.px.back()) = 1.0;
                    seg(sdo.py.back(), sdo.px.back()) = sdo.id.back();
                } else {
                    // Need to assign the source a large segmentation area, let's use half the search radius
                    vec1u idc = where(circular_mask(model.dims,
                        0.5*search_radius/det_aspix, sdo.py.back(), sdo.px.back()) > 0.5);
                    sdo.area.push_back(idc.size());
                    sdo.flux.push_back(total(det_img_clean[idc]));
                    model[idc] = 0;
                    model(sdo.py.back(), sdo.px.back()) = 1.0;
                    for (uint_t i : idc) {
                        if (seg[i] != 0) {
                            uint_t is = where_first(sdo.id == seg[i]);
                            if (is != npos) {
                                sdo.area[is] -= 1.0;
                                sdo.flux[is] -= det_img_clean[i];
                            }
                        }
                    }
                    seg[idc] = sdo.id.back();
                }
            }
        }

        if (opts.verbose) {
            note("found ", sdo.id.size(), " segments in the detection image");
        }

        // Save segmentation map
        if (!opts.nocuts) {
            fits::write(outdir+"segmentation.fits", seg, det_hdr);
            fits::write_table(outdir+"segmentation_catalog.fits", ftable(
                sdo.id, sdo.px, sdo.py, sdo.bx, sdo.by, sdo.area, sdo.origin, sdo.flux
            ));
        }

        // Clean homogenized images
        if (has_homogenize && opts.clean_group) {
            if (opts.verbose) note("cleaning homogenized images...");

            // Make sure the images are homogenized first
            vec2b fitmask = replicate(true, det_img.dims);
            uint_t nhomo = 0;
            for (auto& img : images) {
                if (!img.covered || !img.source.clean_neighbors || img.source.noconvolve) continue;

                if (opts.verbose) note("  preparing ", img.source.filename);

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
            vec2d models(1 + sdo.id.size(), det_img.dims[0]*det_img.dims[1]);
            models(0,_) = 1.0; // background
            vec1u gmodel; gmodel.push_back(0);

            if (opts.clean_model == "psf") {
                // Assume point sources
                for (uint_t s : range(sdo.id)) {
                    models(s+1,_) = flatten(make_point_source_integer(
                        det_psf, det_img.dims, sdo.py[s], sdo.px[s]
                    ));
                }
            } else if (opts.clean_model == "clean") {
                uint_t imod = 1;
                auto c2d = batch_convolve2d(det_psf);

                foreach_segment(seg, sdo.origin, [&](const vec1u& ids) {
                    // Get pixels of that segment
                    vec2d tmod(det_img.dims);
                    tmod[ids] = model[ids]/max(model[ids]);
                    // Convolve to the resolution of the fitted image
                    // NB: we have to assume the same PSF for all homogenized images
                    tmod = c2d.convolve(tmod);
                    // Save model
                    models(imod,_) = flatten(tmod);
                    ++imod;
                });
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

                foreach_segment(seg, sdo.origin, [&](const vec1u& ids) {
                    // Get pixels of that segment
                    vec2d tmod(det_img.dims);
                    tmod[ids] = hdata[ids];
                    double tt = total(tmod[ids])/(rms*sqrt(ids.size()));
                    if (tt > opts.hri_min_snr) {
                        // Keep the HRI model
                        tmod /= max(hdata[ids]);
                        // Convolve to the resolution of the fitted image
                        tmod = c2d.convolve(tmod);
                        models(imod,_) = flatten(tmod);
                    } else {
                        // Not enough SNR, make that a point source
                        models(imod,_) = flatten(make_point_source_integer(
                            det_psf, det_img.dims, sdo.py[imod-1], sdo.px[imod-1]
                        ));
                    }
                    ++imod;
                });
            } else if (opts.clean_model == "nearest_hri") {
                write_error("trying to clean with 'nearest_hri' with group cleaning enabled");
                good = false;
                return;
            }

            // Inspect models and flag out the ones that are zero
            for (uint_t s : range(sdo.id)) {
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
                        write_warning("in the fit for image ", img.source.filename);
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

            if (opts.verbose) note("cleaning ", img.source.filename);

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
            vec2d models(1 + sdo.id.size(), rdata.dims[0]*rdata.dims[1]);
            models(0,_) = 1.0; // background
            vec1u gmodel; gmodel.push_back(0);

            if (opts.clean_model == "psf") {
                // Assume point sources
                for (uint_t s : range(sdo.id)) {
                    models(s+1,_) = flatten(make_point_source_integer(
                        rpsf, rdata.dims, sdo.py[s], sdo.px[s]
                    ));
                }
            } else if (opts.clean_model == "clean") {
                uint_t imod = 1;
                auto c2d = batch_convolve2d(rpsf);

                foreach_segment(seg, sdo.origin, [&](const vec1u& ids) {
                    // Get pixels of that segment
                    vec2d tmod(det_img.dims);
                    tmod[ids] = model[ids]/max(model[ids]);
                    // Convolve to the resolution of the fitted image
                    tmod = c2d.convolve(tmod);
                    // Save model
                    models(imod,_) = flatten(tmod);
                    ++imod;
                });
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

                foreach_segment(seg, sdo.origin, [&](const vec1u& ids) {
                    // Get pixels of that segment
                    vec2d tmod(det_img.dims);
                    tmod[ids] = hdata[ids];
                    double tt = total(tmod[ids])/(rms*sqrt(ids.size()));
                    if (tt > opts.hri_min_snr) {
                        // Keep the HRI model
                        tmod /= max(tmod[ids]);
                        // Convolve to the resolution of the fitted image
                        tmod = c2d.convolve(tmod);
                        models(imod,_) = flatten(tmod);
                    } else {
                        // Not enough SNR, make that a point source
                        models(imod,_) = flatten(make_point_source_integer(
                            rpsf, rdata.dims, sdo.py[imod-1], sdo.px[imod-1]
                        ));
                    }

                    // Save model
                    ++imod;
                });
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
                foreach_segment(seg, sdo.origin, [&](const vec1u& ids) {
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
                            models(imod,_) = flatten(tmod);
                            break;
                        } else {
                            // Not enough SNR, try next HRI
                            ++ihri;
                        }
                    } while (ihri < hri.size());

                    if (ihri == hri.size()) {
                        // Not enough SNR, make that a point source
                        models(imod,_) = flatten(make_point_source_integer(
                            rpsf, rdata.dims, sdo.py[imod-1], sdo.px[imod-1]
                        ));
                    }

                    // Save model
                    ++imod;
                });
            }

            // Inspect models and flag out the ones that are zero
            for (uint_t s : range(sdo.id)) {
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
                    id[im] = sdo.id[gmodel[m]-1];
                    sx[im] = sdo.px[gmodel[m]-1];
                    sy[im] = sdo.py[gmodel[m]-1];

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
}
